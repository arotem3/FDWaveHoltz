#ifndef WAVEHOLTZ_HPP
#define WAVEHOLTZ_HPP

#include <cmath>

#include "Wave.hpp"
#include "ode.hpp"

namespace waveholtz
{
    class WaveHoltz
    {
    private:
        FDWaveEquation wave;
        ode::RungeKutta4<double> rk;
        mutable std::vector<double> _w;

    public:
        WaveHoltz(double dx, int nx, int ny, const char bc[4]) : wave(nx, ny, dx, bc), rk(2*nx*ny), _w(2*nx*ny) {}

        // K(t; omega) := 2/T (cos(omega * t) - 1/4).
        // In the WH iteration we have u^{(n+1)} = \int_0^T K(t) u^{(n)}(t) dt, (T = 2pi/omega)
        static inline double K(double t, double omega)
        {
            return omega / M_PI * (std::cos(omega * t) - 0.25);
        }

        void pi0(double * u, const double omega, const double * F_, double dt=0) const;
        void S(double * u, const double omega, double dt=0) const;
    };

    // selects the number of discretization points needed to resolve the
    // Helmholtz problem with frequency omega one the unit interval [0, 1] using
    // the rule of thumb: h^2 omega^3 == K == constant. By default K = 10.0
    inline int num_points_per_unit_length(double omega, double K = 10.0)
    {
        double h = std::sqrt(K / std::pow(omega, 3));
        return std::ceil(1.0 / h) + 1;
    }
} // namespace waveholtz

#endif