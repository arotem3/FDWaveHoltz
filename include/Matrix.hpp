#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <sstream>

namespace waveholtz
{
    template <typename T>
    class MatrixWrapper
    {
    public:
        const int n_rows;
        const int n_cols;
        T * mat;

        MatrixWrapper(int rows, int cols, T * mat_ = nullptr) : n_rows(rows), n_cols(cols), mat(mat_) {}

        inline const T& operator()(const int row, const int col) const
        {
            return mat[index(row, col)];
        }

        inline T& operator()(const int row, const int col)
        {
            return mat[index(row, col)];
        }

        inline const T& operator()(const int idx) const
        {
            return mat[index(idx)];
        }

        inline T& operator()(const int idx)
        {
            return mat[index(idx)];
        }

    protected:
        inline int index(const int row, const int col) const
        {
            #ifdef WH_DEBUG
            if (row < 0 || row >= n_rows || col < 0 || col >= n_cols)
            {
                std::stringstream err;
                err << "index (" << row << ", " << col << ") out of bounds for matrix of size " << n_rows << " x " << n_cols;
                throw std::out_of_range(err.str());
            }
            #endif
            return row + n_rows * col;
        }

        inline int index(const int idx) const
        {
            #ifdef WH_DEBUG
            if (idx < 0 || idx > n_rows*n_cols)
            {
                std::stringstream err;
                err << "linear index " << idx << " out of bounds for matrix of size " << n_rows << " x " << n_cols << " (= " << n_rows * n_cols << ").";
                throw std::out_of_range(err.str());
            }
            #endif
            return idx;
        }
    };

    template <typename T>
    class Matrix : public MatrixWrapper<T>
    {
    private:
        std::vector<T> _mat;

    public:
        using MatrixWrapper<T>::n_rows;
        using MatrixWrapper<T>::n_cols;
        using MatrixWrapper<T>::mat;

        Matrix(int rows, int cols) : MatrixWrapper<T>(rows, cols), _mat(rows * cols)
        {
            mat = _mat.data();
        }
    };

    typedef Matrix<double> dmat;
    typedef MatrixWrapper<double> dmat_wrapper;
    typedef MatrixWrapper<const double> const_dmat_wrapper;

    template <typename T>
    inline MatrixWrapper<T> reshape(T * data, int nr, int nc)
    {
        return MatrixWrapper<T>(nr, nc, data);
    }
} // namespace waveholtz

#endif