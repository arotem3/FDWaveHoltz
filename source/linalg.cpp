#include "linalg.hpp"

namespace wh
{
    double norm(int n, const double * x)
    {
        double s = 0;

        #pragma omp parallel for reduction(+:s)
        for (int i=0; i < n; ++i)
            s += square(x[i]);
        
        return std::sqrt(s);
    }

    double norm(int n, const std::complex<double> * x)
    {
        double s = 0.0;

        #pragma omp parallel for reduction(+:s)
        for (int i=0; i < n; ++i)
        {
            s += square(x[i].real()) + square(x[i].imag());
        }

        return std::sqrt(s);
    }

    double error(int n, const double * x, const double * y)
    {
        double s = 0;

        #pragma omp parallel for reduction(+:s)
        for (int i=0; i < n; ++i)
            s += square(x[i] - y[i]);
        
        return std::sqrt(s);
    }
} // namespace wh
