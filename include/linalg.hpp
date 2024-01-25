#ifndef LINALG_HPP
#define LINALG_HPP

#include <omp.h>

#include <cmath>
#include <complex>

#include "Tensor.hpp"

namespace wh
{
    inline auto square(auto x)
    {
        return x * x;
    }

    // Euclidean norm of x which has length n.
    double norm(int n, const double * x);
    double norm(int n, const std::complex<double> * x);

    template <int Dim>
    inline double norm(const TensorWrapper<Dim, double>& x)
    {
        return norm(x.size(), x);
    }

    template <int Dim>
    inline double norm(const TensorWrapper<Dim, std::complex<double>>& x)
    {
        return norm(x.size(), x);
    }

    // ||x - y|| in the Euclidean norm where x and y have length n.
    double error(int n, const double * x, const double * y);

    template <int Dim>
    inline double error(const TensorWrapper<Dim, double>& x, const TensorWrapper<Dim, double>& y)
    {
        return error(x.size(), x, y);
    }
} // namespace wh


#endif