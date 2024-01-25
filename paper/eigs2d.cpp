#include <iostream>
#include <iomanip>
#include <fstream>
#include <format>
#include <cmath>

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"
#include "save_binary.hpp"

using namespace wh;

typedef std::complex<double> zdbl;

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

int main()
{
    constexpr int ws[] = {10, 20, 30};

    std::cout << std::setprecision(3);
    std::cout << std::setw(10) << "omega" << " | "
              << std::setw(10) << "#dof" << " | "
              << std::setw(10) << "cond(R)" << " | "
              << std::setw(10) << "time(sec)" << " | "
              << std::endl;

    for (int w : ws)
    {
        Timer stopwatch;

        const double omega = M_PI * w;
        const int nx = num_points_per_unit_length(omega);
        const double h = 2.0 / (nx - 1);

        char boundary_conditions[4];
        boundary_conditions[0] = 'n'; // bottom
        boundary_conditions[1] = 'o'; // right
        boundary_conditions[2] = 'o'; // top
        boundary_conditions[3] = 'n'; // left

        int dims[] = {nx,nx};
        waveholtz2d WH(omega, dims, h, boundary_conditions);
        const int ndof = 2 * nx * nx;

        dmat A(ndof, ndof);
        WH.get_wave_op().as_matrix(A);

        auto [eigvals, R] = eig(A);

        save_binary(eigvals.data(), ndof, std::format("solution/eigs2d_{}", w));

        SingularValueDecomp svd(R);
        std::cout << std::setw(10) << ws << "pi | "
                  << std::setw(10) << ndof << " | "
                  << std::setw(10) << svd.cond() << " | "
                  << std::setw(10) << stopwatch.elapsed() << " | "
                  << std::endl;
    }

    return 0;
}