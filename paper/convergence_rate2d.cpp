#include <iostream>
#include <iomanip>
#include <fstream>
#include <format>
#include <cmath>

#include "WaveHoltz.hpp"
#include "linalg.hpp"
#include "Timer.hpp"

#include "linsolve.hpp"

using namespace wh;

typedef std::complex<double> zdbl;

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
    conv_out << "omega,rho,FP#,GMRES#\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(14) << "ω | "
              << std::setw(13) << "ndof | "
              << std::setw(14) << "ρ | "
              << std::setw(13) << "FP# | "
              << std::setw(13) << "gmres# | "
              << std::setw(13) << "time(sec)"
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

        double err = 1.0;
        bool converged = false; // is ||u[n] - u[n-1]|| / ||pi0|| < tol ?
        int n_iter = max_iter; // # of iterations until ||u[n] - u[n-1]|| / ||pi0|| < tol
        for (int it=2; it <= max_iter; ++it)
        {
            // u_prev = u
            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u_prev[i] = u[i];

            // u = S*u + pi0
            WH.S(u);

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u[i] += pi0[i];

            // err = ||u[it] - u[it-1]|| / ||pi0||
            err = error(ndof, u, u_prev) / pi_zero;

            if (not converged)
            {
                converged = (err < tol);
                if (converged)
                {
                    n_iter = it;
                    break;
                }
            }

            if (save_iters)
            {
                iter_out << err << std::endl;
            }
        }

        const double rho = std::pow(err, 1.0 / n_iter);

        // gmres
        auto IminusS = [&](const double * x, double * y)
        {
            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                y[i] = x[i];

            WH.S(y);

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                y[i] = x[i] - y[i];
        };

        u.zeros();

        linsol::gmres_options<double> options;
        options.absolute_tolerance = 1e-12;
        options.relative_tolerance = tol;
        options.verbose = 0;
        options.restart = 500;
        options.maximum_iterations = 1000;

        auto result = linsol::gmres(ndof, u.data(), IminusS, pi0.data(), options);

        if (save_iters)
        {
            iter_out.close();
            iter_out.open(std::format("solution/gmres2d_{}.txt", (int)w));
            iter_out << std::setprecision(10);
            for (double r : result.residual_norm)
                iter_out << r << std::endl;
            iter_out.close();
        }

        conv_out << omega << ", " << rho << ", " << n_iter << ", " << result.num_matvec << std::endl;
        
        std::cout << std::setw(9) << w << "π | "
                  << std::setw(10) << ndof << " | "
                  << std::setw(10) << rho << " | "
                  << std::setw(10) << n_iter << " | "
                  << std::setw(10) << result.num_matvec << " | "
                  << std::setw(10) << stopwatch.elapsed()
                  << std::endl;
    }

    return 0;
}