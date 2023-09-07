#include "Wave.hpp"

namespace waveholtz
{
    void FDWaveEquation::operator()(double * ut_, const double t, const double * u_) const
    {
        const double one_over_h_squared = 1.0 / h / h;
        const double two_over_h = 2.0 / h;

        auto vt = reshape(ut_ + m*n, m, n);
        auto u = reshape(u_, m, n);
        auto v = reshape(u_+m*n, m, n);

        // du/dt = v
        #pragma omp parallel for
        for (int i=0; i < m*n; ++i)
            ut_[i] = v(i);

        #pragma omp parallel for
        for (int j=1; j < n-1; ++j)
        {
            // BC @ x[0]
            vt(0, j) = one_over_h_squared * (-4.0*u(0, j) + 2.0*u(1, j) + u(0, j-1) + u(0, j+1));

            // interior points
            for (int i=1; i < m-1; ++i)
            {
                vt(i, j) = one_over_h_squared * (-4.0*u(i, j) + u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1));
            }

            // BC @ x[m-1]
            vt(m-1, j) = one_over_h_squared * (-4.0*u(m-1, j) + 2.0*u(m-2, j) + u(m-1, j-1) + u(m-1, j+1));
        }

        // corners
        vt(0, 0) = one_over_h_squared * (-4.0*u(0,0) + 2.0*u(1,0) + 2.0*u(0,1));
        vt(m-1,0) = one_over_h_squared * (-4.0*u(m-1,0) + 2.0*u(m-2,0) + 2.0*u(m-1,1));
        vt(0, n-1) = one_over_h_squared * (-4.0*u(0,n-1) + 2.0*u(1,n-1) + 2.0*u(0,n-2));
        vt(m-1,n-1) = one_over_h_squared * (-4.0*u(m-1,n-1) + 2.0*u(m-2,n-1) + 2.0*u(m-1,n-2));

        #pragma omp parallel for
        for (int i=1; i < m-1; ++i) // BC @ y[0]
            vt(i, 0) = one_over_h_squared * (-4.0*u(i, 0) + u(i-1, 0) + u(i+1, 0) + 2.0*u(i, 1));

        #pragma omp parallel for
        for (int i=1; i < m-1; ++i) // BC @ y[n-1]
            vt(i, n-1) = one_over_h_squared * (-4.0*u(i, n-1) + u(i-1, n-1) + u(i+1, n-1) + 2.0*u(i, n-2));

        if (bc[0] == 'o') // BC @ x[0]
        {
            #pragma omp parallel for
            for (int j=0; j < n; ++j)
                vt(0, j) -= two_over_h * v(0, j);
        }

        if (bc[2] == 'o') // BC @ x[m-1]
        {
            #pragma omp parallel for
            for (int j=0; j < n; ++j)
                vt(m-1, j) -= two_over_h * v(m-1, j);
        }

        if (bc[3] == 'o') // BC @ y[0]
        {
            #pragma omp parallel for
            for (int i=0; i < m; ++i)
                vt(i, 0) -= two_over_h * v(i, 0);
        }

        
        if (bc[1] == 'o') // BC @ y[n-1]
        {
            #pragma omp parallel for
            for (int i=0; i < m; ++i)
                vt(i, n-1) -= two_over_h * v(i, n-1);
        }
    }

} // namespace waveholtz
