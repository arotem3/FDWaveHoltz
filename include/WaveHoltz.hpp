#ifndef WAVEHOLTZ_HPP
#define WAVEHOLTZ_HPP

#include <cmath>
#include <complex>

#include "Wave.hpp"
#include "ode.hpp"

namespace wh
{
    using namespace std::complex_literals;

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
    // the rule of thumb: h^2 omega^3 == K where K is a constant, and K = 10.0 by default.
    inline int num_points_per_unit_length(double omega, double K = 10.0)
    {
        return std::ceil(std::sqrt( std::pow(omega, 3) / K ));
    }

    // returns the scaled filter transfer function beta(z) = 1/pi int_0^{2pi} (cos(t)-1/4)exp(zt) dt
    inline std::complex<double> filter_transfer_function(std::complex<double> z)
    {
        if (z == 0.0)
        {
            return -0.5;
        }
        else if (z == -1.0 || z == 1.0)
        {
            return 0;
        }
        else
        {
            std::complex<double> a = 3.0 * z * z - 1.0;
            std::complex<double> b = std::exp(2 * M_PI * z) - 1.0;
            std::complex<double> c = 4.0 * M_PI * z * (z * z + 1.0);
            return a * b / c;
        }
    }
} // namespace wh

#endif