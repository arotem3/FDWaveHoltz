/* 
    Implements the one dimensional convergence rate analysis.
    
    This file writes:
        solution/convergence_rate.txt
        solution/iters*.txt
    
    convergence_rate.txt is a csv file with header:
        omega,e,mu,min parabolic distance,kappa,avg rate (e),avg rate (mu)
    where omega is the frequency divided by pi and has the values omega_start:omega_delta:omega_end.
    e is ||e_h^(1)|| / ||e_h^(0)||. Similarly, mu is the ||mu_h^(1)|| / ||mu_h^(0)||.
    'min parabolic distance' is the minimum value of -re(lambda/omega) + alpha im(lambda/omega - i) where lambda are the eigenvalues of the discretization matrix.
    kappa is the condition number of the eigenvector matrix.
    'avg rate (e)' is the average value of ||e_h^(n+1)|| / ||e_h^(n)|| for n=1,...,n_iter.
    'avg rate (mu)' is the average value of ||mu_h^(n+1)|| / ||mu_h^(n)|| for n=1,...,n_iter.

    solution/iter*.txt are csv file with * being a value of omega (by default 10, 20, or 30). The header is:
        e,mu
    e is ||e_h^(n)|| / ||e_h^(0)|| for n=1,...,n_iter. Similarly for mu

    Some of this information is printed at runtime to track progress.
*/

#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"

using namespace wh;

typedef std::complex<double> zdbl;

static void initial_condition(double& u, double& v, double x, double w)
{
    const double sp = std::sin(M_PI * x);
    const double s2p = std::sin(2 * M_PI * x);
    const double sw = std::sin(w * x);
    const double cw = std::cos(w * x);
    u =  2.0 * sp * sp * sw;
    v = -2.0 * (w * cw * sp * sp + M_PI * s2p * sw);
}

static dvec linspace(double a, double b, int n)
{
    const double h = (b - a) / (n - 1);
    
    dvec x(n);
    for (int k=0; k < n; ++k)
    {
        x(k) = a + h * k;
    }

    return x;
}

extern "C" void zgesvd_(char * jobu, char * jobvt, int * m, int * n, zdbl * a, int * lda, double * s, zdbl * u, int * ldu, zdbl * vt, int * ldvt, zdbl * work, int * lwork, double * rwork, int * info);
extern "C" void zgemv_(const char *, const int *, const int *, const zdbl *, const zdbl *, const int *, const zdbl *, const int *, const zdbl *, zdbl *, const int *);

class SingularValueDecomp
{
private:
    int n;
    dvec s;
    Matrix<zdbl> U;
    Matrix<zdbl> VH;
    mutable std::vector<zdbl> work;

public:
    // destroys the content of matrix a
    SingularValueDecomp(Matrix<zdbl>& a)
        : n(a.shape()[0]),
          s(n),
          U(n, n),
          VH(n, n),
          work(1)
    {
        char jobu = 'A';
        char jobvt = 'A';
        int lwork = -1;
        std::vector<double> rwork(5 * n);
        int info;

        zgesvd_(&jobu, &jobvt, &n, &n, a, &n, s, U, &n, VH, &n, work.data(), &lwork, rwork.data(), &info);
        lwork = work[0].real();
        work.resize(lwork);

        zgesvd_(&jobu, &jobvt, &n, &n, a, &n, s, U, &n, VH, &n, work.data(), &lwork, rwork.data(), &info);
    }

    double cond() const
    {
        return s(0) / s(n-1);
    }

    // x <- A \ b
    void solve(Vec<zdbl>& x, const auto& b) const
    {
        if (work.size() < n)
            work.resize(n);

        for (int i=0; i < n; ++i)
            work[i] = b[i];
        
        char trans = 'C';
        zdbl zero = 0.0+0.0i;
        zdbl one = 1.0+0.0i;
        int n_ = n;
        int ione = 1;

        zgemv_(&trans, &n, &n, &one, U, &n, work.data(), &ione, &zero, x, &ione);

        for (int i=0; i < n; ++i)
            work[i] = x[i] / s[i];
        
        zgemv_(&trans, &n, &n, &one, VH, &n, work.data(), &ione, &zero, x, &ione);
    }
};

extern "C" void zgeevx_(char * balance, char * jobvl, char * jobvr, char * sense, int * n, zdbl * a, int * lda, zdbl * w, zdbl * vL, int * ldvl, zdbl * vR, int * ldvr, int * ilo, int * ihi, double * scale, double * abnrm, double * rconde, double * rcondv, zdbl * work, int * lwork, double * rwork, int * info);

static std::pair<Vec<zdbl>, Matrix<zdbl>> eig(const dmat& a)
{
    using namespace std::complex_literals;

    int n = a.shape()[0];

    Matrix<zdbl> A(n, n);
    for (int i=0; i < n*n; ++i)
    {
        A[i] = a[i];
    }

    Vec<zdbl> eigenvalues(n);
    Matrix<zdbl> eigenvectors(n, n);

    char balance = 'B';
    char jobvl = 'V';
    char jobvr = 'V';
    char sense = 'B';
    Matrix<zdbl> vl(n, n);
    int ilo, ihi;
    dvec s(n);
    double abnrm;
    dvec rconde(n);
    dvec rcondv(n);
    std::vector<zdbl> work(1);
    dvec rwork(2*n);
    int lwork = -1;
    
    int info;
    // query workspace
    zgeevx_(&balance, &jobvl, &jobvr, &sense, &n, A, &n, eigenvalues, vl, &n, eigenvectors, &n, &ilo, &ihi, s, &abnrm, rconde, rcondv, work.data(), &lwork, rwork, &info);
    // alloc work
    lwork = work[0].real();
    work.resize(lwork);
    // compute eigenvalues and eigenvectors
    zgeevx_(&balance, &jobvl, &jobvr, &sense, &n, A, &n, eigenvalues, vl, &n, eigenvectors, &n, &ilo, &ihi, s, &abnrm, rconde, rcondv, work.data(), &lwork, rwork, &info);

    return make_pair(std::move(eigenvalues), std::move(eigenvectors));
}

// maintain an average value for a sequence
struct running_avg
{
    int n;
    double value;

    running_avg() : n(0), value(0) {}

    // updates the average value and returns it
    double update(double x)
    {
        double a = 1.0 / double(n+1);
        double b = double(n) * a;
        value = a * x + b * value;
        n++;
        return value;
    }
};

int main()
{
    constexpr double alpha = (2.0 * M_PI * M_PI - 3.0) / (12.0 * M_PI);

    const double omega_start = 10, omega_end = 30, omega_delta = 0.5;
    const double a = 0.0, b = 2.0;

    const int n_iter = 1000;

    std::cout << std::fixed << std::setprecision(4);
    
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads." << std::endl;
    }

    std::ofstream out("solution/convergence_rate1d.txt");
    out << std::setprecision(10);
    out << "omega,e,mu,min parabolic distance,kappa,avg rate (e),avg rate (mu),max beta\n";

    std::ofstream iter_out;

    for (double w = omega_start; w <= omega_end; w += omega_delta)
    {
        Timer stopwatch;
        const double omega = w * M_PI;
        const int n = num_points_per_unit_length(omega + 2*M_PI);

        const dvec x = linspace(a, b, n);
        const double h = x(1) - x(0);

        const char bc[] = "no";
        waveholtz1d WH(omega, &n, h, bc);

        dmat A(2*n, 2*n);
        WH.get_wave_op().as_matrix(A);
        auto [eigvals, R] = eig(A);

        double max_beta = 0, min_distance = std::numeric_limits<double>::infinity();
        for (int i=0; i < 2*n; ++i)
        {
            zdbl z = eigvals(i) / omega;

            const double beta = std::abs( filter_transfer_function(z) );
            max_beta = std::max(beta, max_beta);

            const double dist = -z.real() + alpha * square(z.imag() - 1.0);
            min_distance = std::min(dist, min_distance);
        }

        // pre-factor R
        SingularValueDecomp svd(R);
        const double kappa = svd.cond();

        dmat U(n, 2);
        for (int i=0; i < n; ++i)
        {
            initial_condition(U(i, 0), U(i, 1), x(i), omega);
        }

        Vec<zdbl> mu(2*n);
        svd.solve(mu, U);

        const double e0 = norm(U);
        const double mu0 = norm(mu);

        WH.S(U);
        svd.solve(mu, U);

        const double e1 = norm(U);
        const double mu1 = norm(mu);

        const bool save_iters = (int)w % 10 == 0;

        if (save_iters)
        {
            iter_out.open("solution/iters" + std::to_string((int)w) + ".txt");
            iter_out << std::setprecision(10);
            iter_out << "e,mu\n" << e1/e0 << "," << mu1/mu0 << std::endl;
        }

        running_avg r_e, r_mu;
        r_e.update(e1 / e0);
        r_mu.update(mu1 / mu0);

        double e_prev = e1;
        double mu_prev = mu1;
        for (int k = 2; k <= n_iter; ++k)
        {
            WH.S(U);
            const double ek = norm(U);
            
            svd.solve(mu, U);
            const double muk = norm(mu);

            if (save_iters)
            {
                iter_out << ek/e0 << "," << muk/mu0 << std::endl;
            }

            r_e.update(ek / e_prev);
            r_mu.update(muk / mu_prev);

            e_prev = ek;
            mu_prev = muk;
        }

        if (save_iters)
        {
            iter_out.close();
        }

        std::cout << "omega = " << w << " pi | n = " << std::setw(5) << n << " | rate estimate (" << std::setw(8) << 1-min_distance << ") >= max beta (" << std::setw(8) << max_beta << ") >= mu1/mu0 (" << std::setw(8) << mu1 / mu0 << ") | e1 / e0 = " << std::setw(8) << e1 / e0 << " | computation time = " << std::setw(10) << stopwatch.elapsed() << " seconds" << "\n";
        out << w << ", " << e1 / e0 << ", " << mu1 / mu0 << ", " << min_distance << ", " << kappa << ",  " << r_e.value << ", " << r_mu.value << ", " << max_beta << std::endl;
    }

    return 0;
}