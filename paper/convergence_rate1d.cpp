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

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"

#include <armadillo>

using namespace wh;

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

// complex svd factorization for repeated solving linear systems
struct factorization
{
    arma::cx_mat U;
    arma::cx_mat V;
    arma::vec s;

    factorization(const arma::cx_mat& a)
    {
        arma::svd(U, s, V, a);
    }

    arma::cx_vec inv(const arma::vec& x)
    {
        arma::cx_vec y = U.t() * x;
        y = y/s;
        y = V * y;
        return y;
    }

    arma::cx_vec inv(const double* x)
    {
        const int n = U.n_rows;
        arma::vec X(const_cast<double*>(x), n, false, false);
        return inv(X);
    }

    double cond() const
    {
        return s.front() / s.back();
    }
};

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

        arma::mat A(2*n, 2*n);
        WH.get_wave_op().as_matrix( A.memptr() );
        arma::cx_vec eigvals;
        arma::cx_mat R;
        arma::eig_gen(eigvals, R, A, "balance");

        arma::cx_vec beta = eigvals / omega;
        beta.transform(filter_transfer_function);
        const double max_beta = arma::abs(beta).max();

        eigvals /= omega;
        arma::vec parabolic_distance = -arma::real(eigvals) + alpha * arma::square(arma::imag(eigvals) - 1.0);
        const double min_distance = parabolic_distance.min();

        // pre-factor R
        factorization svd(R);
        const double kappa = svd.cond();

        dmat U(n, 2);
        for (int i=0; i < n; ++i)
        {
            initial_condition(U(i, 0), U(i, 1), x(i), omega);
        }

        arma::cx_vec mu = svd.inv(U);

        const double e0 = norm(U);
        const double mu0 = arma::norm(mu);

        WH.S(U);
        mu = svd.inv(U);

        const double e1 = norm(U);
        const double mu1 = arma::norm(mu);

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
            
            mu = svd.inv(U);
            const double muk = arma::norm(mu);

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
        out << w << ", " << e1 / e0 << ", " << mu1 / mu0 << ", " << min_distance << ", " << kappa << ",  " << r_e.value << ", " << r_mu.value << ", " << max_beta << "\n";
    }

    return 0;
}