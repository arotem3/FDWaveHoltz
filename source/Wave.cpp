#include "Wave.hpp"

namespace wh
{
    template <>
    void FDWaveEquation<2>::operator()(double * ut_, const double t, const double * u_) const
    {
        const double one_over_h_squared = 1.0 / h / h;
        const double two_over_h = 2.0 / h;

        const int M = n[0];
        const int N = n[1];

        auto vt = reshape(ut_ + M*N, M, N);
        auto u = reshape(u_, M, N);
        auto v = reshape(u_+M*N, M, N);

        // du/dt = v
        #pragma omp parallel for
        for (int i=0; i < M*N; ++i)
            ut_[i] = v[i];

        #pragma omp parallel for
        for (int j=1; j < N-1; ++j)
        {
            // BC @ x[0]
            vt(0, j) = one_over_h_squared * (-4.0*u(0, j) + 2.0*u(1, j) + u(0, j-1) + u(0, j+1));

            // interior points
            for (int i=1; i < M-1; ++i)
            {
                vt(i, j) = one_over_h_squared * (-4.0*u(i, j) + u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1));
            }

            // BC @ x[m-1]
            vt(M-1, j) = one_over_h_squared * (-4.0*u(M-1, j) + 2.0*u(M-2, j) + u(M-1, j-1) + u(M-1, j+1));
        }

        // corners
        vt(0, 0) = one_over_h_squared * (-4.0*u(0,0) + 2.0*u(1,0) + 2.0*u(0,1));
        vt(M-1,0) = one_over_h_squared * (-4.0*u(M-1,0) + 2.0*u(M-2,0) + 2.0*u(M-1,1));
        vt(0, N-1) = one_over_h_squared * (-4.0*u(0,N-1) + 2.0*u(1,N-1) + 2.0*u(0,N-2));
        vt(M-1,N-1) = one_over_h_squared * (-4.0*u(M-1,N-1) + 2.0*u(M-2,N-1) + 2.0*u(M-1,N-2));

        #pragma omp parallel for
        for (int i=1; i < M-1; ++i) // BC @ y[0]
            vt(i, 0) = one_over_h_squared * (-4.0*u(i, 0) + u(i-1, 0) + u(i+1, 0) + 2.0*u(i, 1));

        #pragma omp parallel for
        for (int i=1; i < M-1; ++i) // BC @ y[n-1]
            vt(i, N-1) = one_over_h_squared * (-4.0*u(i, N-1) + u(i-1, N-1) + u(i+1, N-1) + 2.0*u(i, N-2));

        if (bc[0] == 'o') // BC @ x[0]
        {
            #pragma omp parallel for
            for (int j=0; j < N; ++j)
                vt(0, j) -= two_over_h * v(0, j);
        }

        if (bc[2] == 'o') // BC @ x[m-1]
        {
            #pragma omp parallel for
            for (int j=0; j < N; ++j)
                vt(M-1, j) -= two_over_h * v(M-1, j);
        }

        if (bc[3] == 'o') // BC @ y[0]
        {
            #pragma omp parallel for
            for (int i=0; i < M; ++i)
                vt(i, 0) -= two_over_h * v(i, 0);
        }

        
        if (bc[1] == 'o') // BC @ y[n-1]
        {
            #pragma omp parallel for
            for (int i=0; i < M; ++i)
                vt(i, N-1) -= two_over_h * v(i, N-1);
        }
    }

    template <>
    void FDWaveEquation<2>::as_matrix(double * a) const 
    {
        // lazy!

        int ndof = 2 * _n[0] * _n[1];
        dmat_wrapper A(a, ndof, ndof);

        dvec x(ndof);
        x.zeros();

        const double t = 0.0;

        for (int i=0; i < ndof; ++i)
        {
            x(i) = 1.0;

            operator()(&A(0, i), t, x);

            x(i) = 0.0;
        }
    }

    template <>
    void FDWaveEquation<1>::operator()(double * ut, const double t, const double * u) const
    {
        const double one_over_h_squared = 1.0 / h / h;
        const double two_over_h = 2.0 / h;

        int m = this->n[0];

        double * vt = ut + m;
        const double * v = u + m;

        // du/dt = v
        #pragma omp parallel for
        for (int i=0; i < m; ++i)
            ut[i] = v[i];

        #pragma omp parallel for
        for (int i=1; i < m-1; ++i)
        {
            vt[i] = one_over_h_squared * (u[i-1] - 2.0*u[i] + u[i+1]);
        }
        
        // left BC
        vt[0] = 2.0 * one_over_h_squared * (u[1] - u[0]);
        if (bc[0] == 'o')
        {
            vt[0] -= two_over_h * v[0];
        }

        vt[m-1] = 2.0 * one_over_h_squared * (u[m-2] - u[m-1]);
        if (bc[1] == 'o')
        {
            vt[m-1] -= two_over_h * v[m-1];
        }
    }

    template <>
    void FDWaveEquation<1>::as_matrix(double * a) const
    {
        const double one_over_h_squared = 1.0 / h / h;
        const double two_over_h = 2.0 / h;

        const int m = this->n[0];

        auto A = reshape(a, m, 2, m, 2); // block structure
        
        for (int i=0; i < m; ++i)
            A(i, 0, i, 1) = 1.0;
        
        for (int i=1; i < m-1; ++i)
        {
            A(i, 1, i-1, 0) =  one_over_h_squared;
            A(i, 1, i,   0) = -2.0*one_over_h_squared;
            A(i, 1, i+1, 0) =  one_over_h_squared;
        }

        int i = 0; // left BC
        A(i, 1, i,   0) = -2.0*one_over_h_squared;
        A(i, 1, i+1, 0) =  2.0*one_over_h_squared;
        if (bc[0] == 'o')
            A(i, 1, i, 1) = -two_over_h;

        i = m-1; // right BC
        A(i, 1, i-1, 0) =  2.0*one_over_h_squared;
        A(i, 1, i,   0) = -2.0*one_over_h_squared;
        if (bc[1] == 'o')
            A(i, 1, i, 1) = -two_over_h;
    }
} // namespace wh
