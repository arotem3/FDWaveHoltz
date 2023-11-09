#ifndef WAVE_HPP
#define WAVE_HPP

#include <algorithm>
#include <string>

#include <omp.h>

#include "Tensor.hpp"

namespace wh
{
    template <int Dim>
    class FDWaveEquation
    {
    private:
        const std::string bc;
        int _n[Dim];

    public:
        const double h;
        const int * const n = _n;

        /// @brief initialize wave equation
        /// @param nx number of grid points in x direction
        /// @param ny number of grid points in y direction
        /// @param dx grid spacing, same for x and y
        /// @param bc_ the boundary conditions in the order: bottom, right, top, left
        /// (i.e. clockwise). Either 'n' for Neumann, or 'o' for outflow.

        /// @brief initialize wave equation
        /// @param nx length Dim. Specifies number of grid points along each dimension.
        /// @param dx uniform grid size
        /// @param bc_ length 2*Dim. Specify the boundary conditions. In 1D: bc = {left bc, right bc}. In 2D: bc = {bottom, top, right, left}
        FDWaveEquation(const int nx[], const double dx, const char bc_[]) : bc(bc_), h(dx)
        {
            for (int i=0; i < Dim; ++i)
            {
                _n[i] = nx[i];
            }
        }

        /// @brief evaluates the time derivative of the homogeneous wave equation du/dt = v, dv/dt = laplacian(u).
        /// @param ut (du/dt, dv/dt) has dimension n x 2
        /// @param t time (not used)
        /// @param u (u, v) has dimension n x 2
        void operator()(double * ut, const double t, const double * u) const;
    };

    typedef FDWaveEquation<1> wave1d;
    typedef FDWaveEquation<2> wave2d;
} // namespace wh

#endif