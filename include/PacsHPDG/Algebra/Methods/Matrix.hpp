/**
 * @file Matrix<T>.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MATRIX_METHODS_PACS
#define MATRIX_METHODS_PACS

#include "../Matrix.hpp"

namespace pacs {

    /**
     * @brief LU decomposition of a square matrix.
     * 
     * @tparam T Numeric type.
     * @return std::array<Matrix<T>, 2> 
     */
    template<NumericType T>
    std::array<Matrix<T>, 2> LU(const Matrix<T> &matrix) {
        #ifndef NDEBUG
        assert(matrix.rows == matrix.columns); // Integrity check.
        #endif

        Matrix<T> L{matrix.rows, matrix.columns};
        Matrix<T> U{matrix.rows, matrix.columns};

        for(std::size_t j = 0; j < matrix.columns; ++j) {
            L(j, j) = static_cast<T>(1);

            for(std::size_t i = 0; i <= j; ++i) {
                T sum = static_cast<T>(0);

                for (std::size_t k = 0; k < i; ++k)
                    sum += L(i, k) * U(k, j);

                U(i, j) = matrix(i, j) - sum;
            }

            for(std::size_t i = j + 1; i < matrix.rows; ++i) {
                T sum = static_cast<T>(0);

                for (std::size_t k = 0; k < j; ++k)
                    sum += L(i, k) * U(k, j);

                L(i, j) = (matrix(i, j) - sum) / U(j, j);
            }
        }

        return {L, U};
    }

    /**
     * @brief QR decomposition of a matrix.
     * 
     * @tparam T Numeric type.
     * @return std::array<Matrix<T>, 2> 
     */
    template<NumericType T>
    std::array<Matrix<T>, 2> QR(const Matrix<T> &matrix) {

        Matrix<T> I{matrix.rows, matrix.rows};

        for(std::size_t j = 0; j < matrix.rows; ++j)
            I.elements[j * (matrix.rows + 1)] = static_cast<T>(1);

        Matrix<T> Q{I};
        Matrix<T> R{matrix};

        for(std::size_t j = 0; j < ((matrix.rows > matrix.columns) ? matrix.columns : matrix.rows); ++j) {

            Vector<T> vector{matrix.rows - j};

            for(std::size_t k = 0; k < matrix.rows - j; ++k)
                vector[k] = R.elements[(j + k) * matrix.columns + j];

            vector[0] += (vector[0] > 0 ? static_cast<T>(1) : static_cast<T>(-1)) * norm(vector);
            vector /= norm(vector);

            Matrix<T> H{I};

            for(std::size_t k = 0; k < matrix.rows - j; ++k)
                for(std::size_t l = 0; l < matrix.rows - j; ++l)
                    H.elements[(j + k) * matrix.rows + (j + l)] -= 2 * vector[k] * vector[l];

            R = H * R;
            Q = Q * H.transpose();
        }

        return {Q, R};
    }

    /**
     * @brief Creates an identity matrix.
     * 
     * @tparam T Numeric type.
     * @param size Size of the matrix.
     * @return Matrix<T> 
     */
    template<NumericType T>
    inline Matrix<T> identity(const std::size_t &size) {
        #ifndef NDEBUG // Integrity check.
        assert(size > 0);
        #endif

        Matrix<T> I{size, size};

        for(std::size_t j = 0; j < size; ++j)
            I(j, j) = static_cast<T>(1);

        return I;
    }

    /**
     * @brief Flattens a matrix to a vector.
     * 
     * @tparam T Numeric type.
     * @param matrix Input matrix.
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> squash(const Matrix<T> &matrix) {
        Vector<T> result{matrix.rows * matrix.columns};

        for(std::size_t j = 0; j < matrix.rows; ++j)
            for(std::size_t k = 0; k < matrix.columns; ++k)
                result[j * matrix.columns + k] = matrix(j, k);

        return result;
    }

    /**
     * @brief Product of diagonal elements.
     * 
     * @tparam T Numeric type.
     * @param matrix Input matrix.
     * @return T 
     */
    template<NumericType T>
    T mtrace(const Matrix<T> &matrix) {
        #ifndef NDEBUG // Integrity check.
        assert(matrix.rows == matrix.columns);
        #endif

        T product = static_cast<T>(1);

        for(std::size_t j = 0; j < matrix.rows; ++j)
            product *= matrix(j, j);

        return product;
    }
    
}

#endif