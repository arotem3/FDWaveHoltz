#ifndef LINALG_HPP
#define LINALG_HPP

#include <omp.h>

#include <cmath>

namespace wh
{
    inline double square(double x)
    {
        return x * x;
    }

    // Euclidean norm of x which has length n.
    double norm(int n, const double * x);

    // ||x - y|| in the Euclidean norm where x and y have length n.
    double error(int n, const double * x, const double * y);
} // namespace wh


#endif