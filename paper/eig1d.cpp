#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"

using namespace wh;

typedef std::complex<double> zdbl;

extern "C" void zgeevx_(char * balance, char * jobvl, char * jobvr, char * sense, int * n, zdbl * a, int * lda, zdbl * w, zdbl * vL, int * ldvl, zdbl * vR, int * ldvr, int * ilo, int * ihi, double * scale, double * abnrm, double * rconde, double * rcondv, zdbl * work, int * lwork, double * rwork, int * info);

static Vec<zdbl> eig(const dmat& a)
{
    using namespace std::complex_literals;

    int n = a.shape()[0];

    Matrix<zdbl> A(n, n);
    for (int i=0; i < n*n; ++i)
    {
        A[i] = a[i];
    }

    Vec<zdbl> eigenvalues(n);
    // Matrix<zdbl> eigenvectors(n, n);

    char balance = 'B';
    char jobvl = 'N';
    char jobvr = 'N';
    char sense = 'N';
    // Matrix<zdbl> vl(n, n);
    zdbl * vl = nullptr, * vr = nullptr;
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
    zgeevx_(&balance, &jobvl, &jobvr, &sense, &n, A, &n, eigenvalues, vl, &n, vr, &n, &ilo, &ihi, s, &abnrm, rconde, rcondv, work.data(), &lwork, rwork, &info);
    // alloc work
    lwork = work[0].real();
    work.resize(lwork);
    // compute eigenvalues and eigenvectors
    zgeevx_(&balance, &jobvl, &jobvr, &sense, &n, A, &n, eigenvalues, vl, &n, vr, &n, &ilo, &ihi, s, &abnrm, rconde, rcondv, work.data(), &lwork, rwork, &info);

    return eigenvalues;
}

int main()
{
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads." << std::endl;
    }

    constexpr double C = 10.0;

    const double a = 0.0, b = 2.0;

    for (int w : {10, 20, 30})
    {
        double omega = w * M_PI;

        const int n = (b - a) * num_points_per_unit_length(omega + 2*M_PI, C);
        const double h = (b - a) / (n - 1);
        const char bc[] = "no";

        std::cout << w << "pi | n = " << n << " | h = " << h << " | #dof = " << 2*n << "\n";

        wave1d op(&n, h, bc);
        dmat A(2*n, 2*n);
        op.as_matrix(A);

        auto eigvals = eig(A);

        std::ofstream out(std::format("solution/eig_{}.txt", w));
        out << std::setprecision(10) << std::scientific;
        for (int i=0; i < 2*n; ++i)
        {
            out << std::setw(25) << eigvals.at(i).real() << std::setw(25) << eigvals.at(i).imag() << "\n";
        }
        out.close();
    }

    return 0;
}