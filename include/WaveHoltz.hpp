#ifndef WAVEHOLTZ_HPP
#define WAVEHOLTZ_HPP

#include <cmath>

#include "Wave.hpp"
#include "ode.hpp"

namespace wh
{
    template <int Dim>
    class WaveHoltz
    {
    private:
        int ndof;
        double T;
        double dt;
        int nt;
        FDWaveEquation<Dim> wave;
        ode::RungeKutta4<double> rk;
        mutable dvec w;

    public:
        const double omega;

        WaveHoltz(double omega_, const int nx[], double dx, const char bc[], double dt=0.0);

        // K(t; omega) := 2/T (cos(omega * t) - 1/4).
        // In the WH iteration we have u^{(n+1)} = \int_0^T K(t) u^{(n)}(t) dt, (T = 2pi/omega)
        inline double K(double t) const
        {
            return omega / M_PI * (std::cos(omega * t) - 0.25);
        }

        void pi0(double * u, const double * F) const;
        
        void S(double * u) const;
        
        // map (u, v) to (Re{U}, Im{U}) where U is the solution to Helmholtz if (u, v) = S(u, v) + pi0.
        void postprocess(double * u) const
        {
            const int nn = ndof / 2;

            double * v = u + nn;

            #pragma omp parallel for
            for (int i=0; i < nn; ++i)
            {
                v[i] /= -omega;
            }
        }
    
        const FDWaveEquation<Dim>& get_wave_op() const
        {
            return wave;
        }
    };

    typedef WaveHoltz<1> waveholtz1d;
    typedef WaveHoltz<2> waveholtz2d;

    // selects the number of discretization points needed to resolve the
    // Helmholtz problem with frequency omega one the unit interval [0, 1] using
    // the rule of thumb: h^2 omega^3 == K == constant. By default K = 10.0
    inline int num_points_per_unit_length(double omega, double K = 10.0)
    {
        double h = std::sqrt(K / std::pow(omega, 3));
        return std::ceil(1.0 / h) + 1;
    }
} // namespace wh

#endif