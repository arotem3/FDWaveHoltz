#include "WaveHoltz.hpp"

namespace waveholtz
{
    void WaveHoltz::pi0(double * u, const double omega, const double * F_, double dt) const
    {
        const int nx = wave.m;
        const int ny = wave.n;
        const int ndof = 2 * nx * ny;

        auto F = reshape(F_, nx, ny);

        auto time_derivative = [&](double * p, const double t, const double * y) -> void
        {
            wave(p, t, y); // du/dt = v, dv/dt = laplacian(u)

            // add forcing
            const double s = std::cos(omega * t); // std::exp(1i * omega * t);
            double * vt = p + nx*ny;

            #pragma omp parallel for
            for (int i=0; i < nx*ny; ++i)
                vt[i] += s * F(i);
        };

        const double T = 2.0 * M_PI / omega;
        if (dt <= 0)
            dt = 0.1 * wave.h;
        const int nt = std::ceil(T / dt);
        dt = T / nt;

        double t = 0.0;

        double * w = _w.data();
        
        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
        {
            u[i] = 0.0;
            w[i] = 0.0;
        }

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, time_derivative, t, w);

            double Kh = K(t, omega) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u[i] += Kh * w[i];
        }
    }

    void WaveHoltz::S(double * u, const double omega, double dt) const
    {
        const int ndof = 2 * wave.m * wave.n;
        const double T = 2.0 * M_PI / omega;
        if (dt <= 0.0)
            dt = 0.1 * wave.h;
        const int nt = std::ceil(T / dt);
        dt = T / nt;

        double t = 0.0;
        double * w = _w.data();

        double Kh = 0.5 * K(t, omega) * dt;
        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w[i] = Kh * u[i];

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, wave, t, u);

            double Kh = K(t, omega) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                w[i] += Kh * u[i];
        }

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            u[i] = w[i];
    }
} // namespace waveholtz