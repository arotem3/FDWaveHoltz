#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include <omp.h>

#include "Wave.hpp"
#include "ode.hpp"
#include "Timer.hpp"
#include "save_binary.hpp"
#include "ProgressBar.hpp"

using namespace wh;

static inline double f(double x, double y)
{
    double r = x*x + y*y;
    return 100 * std::exp(-100*r);
}

static inline double initial_value(double x, double y)
{
    return 0.0;
}

static inline double initial_velocity(double x, double y)
{
    return 0.0;
}

int main()
{
    constexpr int save_solution = 20; // save solution to file every save_solution time steps. If save_solution == 0, then solutions will not be saved to file.
    const char solution_base_file_name[] = "solution/u"; // solutions will be of form: solution_base_file_name00i where i is the iteration padded with zeros, e.g. "solution/u0021"
    constexpr int zero_pad = 6;

    auto save = [=](MatrixWrapper<double>& x, int it) -> void
    {
        if (save_solution && (it % save_solution == 0))
        {
            std::stringstream fname;
            fname << solution_base_file_name << std::setw(zero_pad) << std::setfill('0') << it;

            save_binary(x, x.size(), fname.str());
        }
    };

    constexpr int n = 250, ndof = 2 * n * n;
    constexpr double a = -1., b = 1.;
    constexpr double h = (b - a) / (n - 1);

    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads\n";
    }

    dvec w(ndof);
    MatrixWrapper<double> u(w, n, n);
    MatrixWrapper<double> v(w + n*n, n, n);

    dvec x(n); // grid points
    #pragma omp parallel for
    for (int i=0; i < n; ++i)
        x[i] = a + i*h;

    #pragma omp parallel for
    for (int i=0; i < n; ++i)
    {
        for (int j=0; j < n; ++j)
        {
            u(i, j) = initial_value(x[i], x[j]);
            v(i, j) = initial_velocity(x[i], x[j]);
        }
    }

    save(u, 0); // save solution

    dmat F(n, n);
    #pragma omp parallel for
    for (int j=0; j < n; ++j)
        for (int i=0; i < n; ++i)
            F(i, j) = f(x[i], x[j]);

    char boundary_conditions[4];
    boundary_conditions[0] = 'n'; // bottom
    boundary_conditions[1] = 'o'; // right
    boundary_conditions[2] = 'o'; // top
    boundary_conditions[3] = 'n'; // left

    int dims[] = {n,n};
    wave2d wave(dims, h, boundary_conditions);

    auto time_derivative = [&](double * p, double t, const double * y) -> void
    {
        wave(p, t, y); // du/dt = v, dv/dt = laplacian(u)

        // add forcing
        const double s = std::sin(30*t);
        double * vt = p + n*n;

        #pragma omp parallel for
        for (int i=0; i < n*n; ++i)
            vt[i] += s * F[i];
    };

    double dt = 2.5e-4;
    const double T = 2.0;
    const int nt = std::ceil(T / dt);
    dt = T / nt;

    double t = 0.0;
    ode::RungeKutta2<double> rk(ndof);

    std::cout << "n := " << n << "\n"
              << "ndof := " << ndof << "\n"
              << "h := " << h << "\n"
              << "T := " << T << "\n"
              << "dt := " << dt << "\n"
              << "nt := " << nt << "\n";

    Timer stopwatch;
    ProgressBar progress_bar(nt);

    for (int it=1; it <= nt; ++it)
    {
        rk.step(dt, time_derivative, t, w);

        save(u, it);

        ++progress_bar;
        std::cout << "[" << progress_bar.get_progress() << "]" << std::setw(5) << it << " / " << nt << "\r" << std::flush;
    }

    std::cout << "\nTotal time: " << stopwatch.elapsed() << " seconds\n";

    return 0;
}