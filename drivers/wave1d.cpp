/* 
    In this file we solve:
        u_{tt} = u_{xx} + f(x, t)

    with boundary conditions either:
        Neumann: n.u_x = 0
        Outflow: u_t + n.u_x = 0
    
        where n is the unit normal at the boundary.
        So n.u_x is the directional derivative in the direction of n.
 */


#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "Wave.hpp"
#include "ode.hpp"
#include "Timer.hpp"
#include "save_binary.hpp"
#include "ProgressBar.hpp"

using namespace waveholtz;

static inline double f(double x, double t)
{
    return 100 * std::exp(-100*x*x) * std::cos(30 * t);
}

static inline double initial_displacement(double x)
{
    return 0.0;
}

static inline double initial_velocity(double x)
{
    return 0.0;
}

int main()
{
    constexpr int save_solution = 10; // save solution to file every save_solution time steps.
                                      // If save_solution == 0, then solutions will not be saved to file.
    const char solution_base_file_name[] = "solution/u"; // solutions will be of form: solution_base_file_name00i
                                                         // where i is the iteration padded with zeros.
                                                         // e.g. "solution/u0021"
    constexpr int zero_pad = 6;

    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "using " << omp_get_num_threads() << " threads\n";
    }

    auto save = [=](int n, const double * x, int it) -> void
    {
        if (save_solution && it % save_solution == 0)
        {
            std::stringstream fname;
            fname << solution_base_file_name << std::setw(zero_pad) << std::setfill('0') << it;

            save_binary(x, n, fname.str());
        }
    };

    constexpr int n = 250, ndof = 2 * n;
    constexpr double a = -1., b = 1.;
    constexpr double h = (b - a) / (n - 1);

    std::vector<double> w(ndof, 0.0);
    double * u = w.data();
    double * v = u + n;

    std::vector<double> x(n);
    
    #pragma omp parallel for
    for (int i=0; i < n; ++i)
    {
        x[i] = a + i*h;
        u[i] = initial_displacement(x[i]);
        v[i] = initial_velocity(x[i]);
    }

    save_binary(x.data(), n, "solution/x");
    save(n, u, 0);

    char boundary_conditions[2];
    boundary_conditions[0] = 'n'; // bc @ x == a
    boundary_conditions[1] = 'o'; // bc @ x == b

    FDWaveEquation1D wave(n, h, boundary_conditions);

    auto time_derivative = [&](double * p, double t, const double * y) -> void
    {
        wave(p, t, y);

        double * vt = p + n;

        #pragma omp parallel for
        for (int i=0; i < n; ++i)
        {
            vt[i] += f(x[i], t);
        }
    };

    double dt = 2.5e-4;
    double T = 2.0;
    const int nt = std::ceil(T/dt);
    dt = T / nt;

    double t = 0.0;
    ode::RungeKutta2<double> rk(ndof);

    std::cout << "n = " << n << "\n"
              << "h = " << h << "\n"
              << "dt = " << dt << "\n"
              << "nt = " << nt << "\n";

    Timer stopwatch;
    ProgressBar progress_bar(nt);

    for (int it=1; it <= nt; ++it)
    {
        rk.step(dt, time_derivative, t, u);
        
        save(n, u, it);

        ++progress_bar;
        std::cout << "[" << progress_bar.get_progress() << "]" << std::setw(5) << it << " / " << nt << "\r" << std::flush; 
    }

    std::cout << "\nTotal time: " << stopwatch.elapsed() << " seconds\n";

    return 0;
}