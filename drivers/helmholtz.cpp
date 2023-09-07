#include "WaveHoltz.hpp"
#include "Timer.hpp"
#include "ProgressBar.hpp"
#include "save_binary.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

#include <omp.h>

using namespace waveholtz;

inline double square(double x)
{
    return x * x;
}

double error(int n, double * x, double * y)
{
    double s = 0;

    #pragma omp parallel for reduction(+:s)
    for (int i=0; i < n; ++i)
        s += square(x[i] - y[i]);
    
    return std::sqrt(s);
}

int main()
{
    constexpr double K = 20.0;
    constexpr double omega = 30 * M_PI;
    constexpr double tol = 1e-5;
    constexpr int maxit = 1'000;

    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads\n";
    }

    auto f = [omega](double x, double y) -> double
    {
        double r = square(omega * (x + 0.7)) + square(omega * (y+0.1));
        return square(omega) * std::exp(-r);
    };

    const int n = 2 * num_points_per_unit_length(omega, K);
    const double h = 2.0 / (n-1);
    const double dt = 0.1 * h;

    std::cout << "n := " << n << "\n";

    // specify grid [-1, 1] x [-1, 1]
    std::vector<double> x(n);
    const std::vector<double>& y = x;

    #pragma omp parallel for
    for (int i=0; i < n; ++i)
        x[i] = -1.0 * h*i;

    // compute source function F
    dmat F(n, n);
    #pragma omp parallel for collapse(2)
    for (int j=0; j < n; ++j)
        for (int i=0; i < n; ++i)
            F(i, j) = f(x[i], y[j]);

    // specify initial conditions
    char boundary_conditions[4];
    boundary_conditions[0] = 'n'; // bottom
    boundary_conditions[1] = 'o'; // right
    boundary_conditions[2] = 'o'; // top
    boundary_conditions[3] = 'n'; // left

    // WaveHoltz iteration operator
    WaveHoltz Pi(h, n, n, boundary_conditions);
    const int ndof = 2 * n * n;

    std::vector<double> w(ndof, 0.0);
    std::vector<double> w_prev(ndof);
    std::vector<double> pi0(ndof);

    Timer stopwatch;

    // pi0 = \int_0^T K(t) q(t) dt
    // where:
    // q_t = A q - F(x, y) cos(omega t),
    // q(0) = 0,
    // q_t(0) = 0.
    Pi.pi0(pi0.data(), omega, F.mat, dt);

    double pi_zero = 0; // ||pi_zero||
    #pragma omp parallel for reduction(+:pi_zero)
    for (int i=0; i < ndof; ++i)
        pi_zero += square(pi0[i]);
    pi_zero = std::sqrt(pi_zero);

    std::cout << std::scientific << std::setprecision(2);

    // iteration
    for (int it=1; it <= maxit; ++it)
    {
        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w_prev[i] = w[i];

        // S*w = \int_0^T K(t) q(t) dt,
        // where:
        // q_t = A q
        // q(0) = w
        Pi.S(w.data(), omega, dt);

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w[i] += pi0[i];

        double err = error(ndof, w.data(), w_prev.data());

        std::cout << std::setw(10) << it << " / " << maxit << " || error =" << std::setw(10) << err / pi_zero << "\r" << std::flush;

        if (err < tol * pi_zero)
            break;
    }
    std::cout << "\nComputation time: " << stopwatch.elapsed() << " seconds\n";

    save_binary(w.data(), 2*n*n, "waveholtz_solution.dat");

    return 0;
}