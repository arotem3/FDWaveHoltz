#include <iostream>
#include <iomanip>
#include <fstream>
#include <format>
#include <cmath>

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"

using namespace wh;

typedef std::complex<double> zdbl;

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

static double force(double x, double y, double omega)
{
    const double r = square(x + 0.7) + square(y + 0.1);
    const double s = square(omega);
    return (s / M_PI) * std::exp(-0.5 * s * r);
}

int main()
{
    const int max_iter = 1'000;
    const double tol = 1e-6;

    const double omega_start = 10, omega_end = 30, omega_delta = 0.5;

    std::ofstream conv_out("solution/convergence_rate2d.txt");
    conv_out << "omega,r1,avg(r),iter\n";

    std::cout << std::setprecision(3);
    std::cout << std::setw(10) << "omega" << " | "
              << std::setw(10) << "ndof" << " | "
              << std::setw(10) << "r[1]" << " | "
              << std::setw(10) << "avg(r)" << " | "
              << std::setw(10) << "#iter" << " | "
              << std::setw(10) << "time(sec)"
              << std::endl;

    for (double w = omega_start; w <= omega_end; w += omega_delta)
    {
        Timer stopwatch;

        const double omega = M_PI * w;
        const int nx = num_points_per_unit_length(omega);
        const double h = 2.0 / (nx - 1);

        dvec x(nx);
        #pragma omp parallel for
        for (int i=0; i < nx; ++i)
            x(i) = -1.0 + h*i;
        
        dmat F(nx, nx);
        #pragma omp parallel for collapse(2)
        for (int j=0; j < nx; ++j)
            for (int i=0; i < nx; ++i)
                F(i, j) = force(x(i), x(j), omega);
        
        char boundary_conditions[4];
        boundary_conditions[0] = 'n'; // bottom
        boundary_conditions[1] = 'o'; // right
        boundary_conditions[2] = 'o'; // top
        boundary_conditions[3] = 'n'; // left

        int dims[] = {nx,nx};
        waveholtz2d WH(omega, dims, h, boundary_conditions);
        const int ndof = 2 * nx * nx;

        dvec u(ndof);
        u.zeros();
        
        dvec u_prev(ndof);

        dvec pi0(ndof);

        WH.pi0(pi0, F);
        const double pi_zero = norm(ndof, pi0);

        const bool save_iters = (int(2*w) % 20 == 0); // 10, 20, or 30
        std::ofstream iter_out;
        if (save_iters)
        {
            iter_out.open(std::format("solution/iter2d_{}.txt", (int)w));
        }

        double err_prev, err = 1.0;
        running_avg ratio;
        double r1; // ||u[2] - u[1]|| / ||pi0|| = ||u[2] - pi0|| / ||pi0||
        bool converged = false; // is ||u[n] - u[n-1]|| / ||pi0|| < tol ?
        int n_iter = max_iter; // # of iterations until ||u[n] - u[n-1]|| / ||pi0|| < tol
        for (int it=1; it <= max_iter; ++it)
        {
            // u_prev = u
            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u_prev[i] = u[i];

            err_prev = err;

            // u = S*u + pi0
            WH.S(u);

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u[i] += pi0[i];

            // err = ||u[it] - u[it-1]|| / ||pi0||
            err = error(ndof, u, u_prev) / pi_zero;

            // r = ||u[it] - u[it-1]|| / ||u[it-1] - u[it-2]||
            ratio.update(err / err_prev);

            if (it == 2)
                r1 = err;

            if (not converged)
            {
                converged = (err < tol);
                if (converged)
                    n_iter = it;
            }

            if (save_iters)
            {
                iter_out << err << std::endl;
            }
        }

        conv_out << omega << ", " << r1 << ", " << ratio.value << ", " <<  n_iter << std::endl;
        
        std::cout << std::setw(10) << omega << " | "
                  << std::setw(10) << ndof << " | "
                  << std::setw(10) << r1 << " | "
                  << std::setw(10) << ratio.value << " | "
                  << std::setw(10) << n_iter << " | "
                  << std::setw(10) << stopwatch.elapsed()
                  << std::endl;
    }

    return 0;
}