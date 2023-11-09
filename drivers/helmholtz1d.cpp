#include "WaveHoltz.hpp"
#include "Timer.hpp"
#include "ProgressBar.hpp"
#include "save_binary.hpp"
#include "linalg.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace wh;

int main()
{
    constexpr double K = 1.0;
    constexpr double omega = 20;
    constexpr double tol = 1e-8;
    constexpr int maxit = 10'000;

    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads\n";
    }

    auto f = [omega](double x) -> double
    {
        double r = square(omega * (x + 0.7));
        return square(omega) * std::exp(-r);
    };

    const int n = 2 * num_points_per_unit_length(omega, K);
    const int ndof = 2 * n;
    const double h = 2.0 / (n-1);

    std::cout << "n := " << n << "\n";

    dvec x(n); // grid [-1, 1]
    dvec F(n); // source function F

    #pragma omp parallel for
    for (int i=0; i < n; ++i)
    {
        x[i] = -1.0 + h*i;
        F[i] = f(x[i]);
    }
    
    // specify boundary conditions
    char boundary_conditions[2];
    boundary_conditions[0] = 'n'; // left
    boundary_conditions[1] = 'o'; // right

    // WaveHoltz iteration operator
    waveholtz1d WH(omega, &n, h, boundary_conditions);

    dvec w(ndof);
    dvec w_prev(ndof);
    dvec pi0(ndof);

    Timer stopwatch;

    // pi0 = \int_0^T K(t) q(t) dt
    // where:
    // dq/dt = A q - F(x) cos(omega t),
    //  q(0) = 0,
    //  q_t(0) = 0.
    WH.pi0(pi0, F);

    double pi_zero = norm(ndof, pi0); // ||pi_zero||

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
        WH.S(w);

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w[i] += pi0[i];

        double err = error(ndof, w, w_prev);

        std::cout << std::setw(10) << it << " / " << maxit << " || error =" << std::setw(10) << err / pi_zero << "\r" << std::flush;

        if (err < tol * pi_zero)
            break;
    }
    WH.postprocess(w);
    std::cout << "\nComputation time: " << stopwatch.elapsed() << " seconds\n";

    save_binary(x, n, "solution/x");
    save_binary(w, ndof, "solution/u");

    return 0;
}