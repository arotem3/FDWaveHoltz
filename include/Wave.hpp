#ifndef WAVE_HPP
#define WAVE_HPP

#include <algorithm>
#include <string>

#include <omp.h>

#include "Matrix.hpp"

namespace waveholtz
{
    class FDWaveEquation
    {
    private:
        const std::string bc;

    public:
        const double h;
        const int m;
        const int n;

        /// @brief initialize wave equation
        /// @param nx number of grid points in x direction
        /// @param ny number of grid points in y direction
        /// @param dx grid spacing, same for x and y
        /// @param bc_ the boundary conditions in the order: bottom, right, top, left
        /// (i.e. clockwise). Either 'n' for Neumann, or 'o' for outflow. 
        FDWaveEquation(const int nx, const int ny, const double dx, const char bc_[4]) : bc(bc_), h(dx), m(nx), n(ny) {}

        /// @brief evaluates the time derivative of the homogeneous wave equation du/dt = v, dv/dt = laplacian(u).
        /// @param ut (du/dt, dv/dt) has dimension m x n x 2
        /// @param t time (not used)
        /// @param u (u, v) has dimension m x n x 2
        void operator()(double * ut, const double t, const double * u) const;
    };
} // namespace waveholtz

#endif