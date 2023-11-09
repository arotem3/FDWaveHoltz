#include "WaveHoltz.hpp"

namespace wh
{
    template <>
    WaveHoltz<2>::WaveHoltz(double omega_, const int nx[], double dx, const char bc[], double dt_) : ndof(2 * nx[0] * nx[1]), wave(nx, dx, bc), rk(ndof), w(ndof), omega(omega_)
    {
        T = 2 * M_PI / omega;
        dt = (dt_ > 0) ? dt_ : (0.1 * dx);
        nt = std::ceil(T/dt);
        dt = T / nt;
    }

    template <>
    void WaveHoltz<2>::pi0(double * u, const double * F_) const
    {
        const int nx = wave.n[0];
        const int ny = wave.n[1];

        auto F = reshape(F_, nx, ny);

        auto time_derivative = [&](double * p, const double t, const double * y) -> void
        {
            wave(p, t, y); // du/dt = v, dv/dt = laplacian(u)

            // add forcing
            const double s = std::cos(omega * t); // std::exp(1i * omega * t);
            double * vt = p + nx*ny;

            #pragma omp parallel for
            for (int i=0; i < nx*ny; ++i)
                vt[i] -= s * F[i];
        };

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
        {
            u[i] = 0.0;
            w[i] = 0.0;
        }

        double t = 0.0;

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, time_derivative, t, w);

            double Kh = K(t) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u[i] += Kh * w[i];
        }
    }

    template <>
    void WaveHoltz<2>::S(double * u) const
    {
        double t = 0.0;

        double Kh = 0.5 * K(t) * dt;
        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w[i] = Kh * u[i];

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, wave, t, u);

            double Kh = K(t) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                w[i] += Kh * u[i];
        }

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            u[i] = w[i];
    }

    template <>
    WaveHoltz<1>::WaveHoltz(double omega_, const int nx[], double dx, const char bc[], double dt_) : ndof(2 * nx[0]), wave(nx, dx, bc), rk(ndof), w(ndof), omega(omega_)
    {
        T = 2 * M_PI / omega;
        dt = (dt_ > 0) ? dt_ : (0.1 * dx);
        nt = std::ceil(T/dt);
        dt = T / nt;
    }

    template <>
    void WaveHoltz<1>::pi0(double * u, const double * F) const
    {
        const int n = wave.n[0];

        auto time_derivative = [&](double * p, const double t, const double * y) -> void
        {
            wave(p, t, y); // du/dt = v, dv/dt = laplacian(u)

            // add forcing
            const double s = std::cos(omega * t); // std::exp(1i * omega * t);
            double * vt = p + n;

            #pragma omp parallel for
            for (int i=0; i < n; ++i)
                vt[i] -= s * F[i];
        };

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
        {
            u[i] = 0.0;
            w[i] = 0.0;
        }

        double t = 0.0;

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, time_derivative, t, w);

            double Kh = K(t) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                u[i] += Kh * w[i];
        }
    }

    template <>
    void WaveHoltz<1>::S(double * u) const
    {
        double t = 0.0;

        double Kh = 0.5 * K(t) * dt;
        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            w[i] = Kh * u[i];

        for (int it=1; it <= nt; ++it)
        {
            rk.step(dt, wave, t, u);

            double Kh = K(t) * dt;
            Kh *= (it == nt) ? 0.5 : 1.0;

            #pragma omp parallel for
            for (int i=0; i < ndof; ++i)
                w[i] += Kh * u[i];
        }

        #pragma omp parallel for
        for (int i=0; i < ndof; ++i)
            u[i] = w[i];
    }
} // namespace wh