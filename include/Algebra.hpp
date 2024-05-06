/**
 * @file Algebra.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_PACS
#define ALGEBRA_PACS

// Type.
#include <Type.hpp>

// Vectors.
#include <Vector.hpp>

// Matrices.
#include <Matrix.hpp>
#include <Sparse.hpp>

// Concepts.
#include <concepts>

// OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// Algebra tolerance.
#ifndef ALGEBRA_TOLERANCE
#define ALGEBRA_TOLERANCE 1E-8
#endif

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 512
#endif

namespace pacs {

    // METHODS.

    /**
     * @brief Vector Euclidean norm.
     * 
     * @tparam T 
     * @param vector 
     * @return double 
     */
    template<NumericType T>
    double norm(const Vector<T> &vector) {
        double sum = 0.0;

        #pragma omp parallel for reduction(+: sum)
        for(std::size_t j = 0; j < vector.length; ++j)
            sum += std::abs(vector[j]) * std::abs(vector[j]);

        return std::sqrt(sum);
    }

    /**
     * @brief Linear system solver using CGM.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @return Vector<T> 
     */
    template<MatrixLike M, NumericType T>
    Vector<T> solve(const M &A, const Vector<T> &b) {
        #ifndef NDEBUG // Integrity check.
        assert(A.columns == b.length);
        #endif

        Vector<T> xk{b.length};
        Vector<T> rk = b, rkk = b, pk = rk;
        T ak, bk;

        std::size_t iterations = 0;

        // *k: k-th value.
        // *kk: (k + 1)-th value
        do {
            // Residual.
            rk = rkk;

            // New guess point.
            ak = (rk * rk) / (pk * (A * pk));
            xk += ak * pk;

            // New residual.
            rkk -= ak * A * pk;

            // bk and pk.
            bk = (rkk * rkk) / (rk * rk);
            pk = rkk + bk * pk;

            ++iterations;
        } while((norm(rkk) > ALGEBRA_TOLERANCE) && (iterations < ALGEBRA_ITER_MAX));

        return xk;
    }
}

#endif