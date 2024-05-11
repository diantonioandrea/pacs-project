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
#define ALGEBRA_ITER_MAX 2E4
#endif

namespace pacs {

    // METHODS.

    /**
     * @brief Vector Euclidean norm.
     * 
     * @tparam T 
     * @param vector 
     * @return Real 
     */
    template<NumericType T>
    Real norm(const Vector<T> &vector) {
        Real sum = 0.0;

        #pragma omp parallel for reduction(+: sum)
        for(std::size_t j = 0; j < vector.length; ++j)
            sum += std::abs(vector[j]) * std::abs(vector[j]);

        return std::sqrt(sum);
    }

    /**
     * @brief Linear system solver using Gauss-Seidel.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @return Vector<T> 
     */
    template<MatrixLike M, NumericType T>
    Vector<T> solve(const M &matrix, const Vector<T> &vector) {
        #ifndef NDEBUG // Integrity checks.
        assert(matrix.rows == matrix.columns);
        assert(matrix.columns == vector.length);
        #endif

        #ifdef VERBOSE
        std::cout << "Solving a linear system." << std::endl;
        #endif

        // Matrix size.
        const std::size_t size = matrix.rows;

        // Iterations.
        std::size_t iterations = 0;

        // L, diagonal and U.
        M lower = matrix.lower() + matrix.diagonal();
        M upper = matrix.upper();

        // (L + diagonal)'s determinant.
        T determinant = mtrace(lower);
        
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(determinant) > ALGEBRA_TOLERANCE);
        #endif

        // L's inverse.
        M lower_inv{size, size};

        for(std::size_t j = 0; j < size; ++j) {
            Vector<T> column{size};

            for(std::size_t k = 0; k < size; ++k) {
                T sum = static_cast<T>(0);

                for(std::size_t l = 0; l < k; ++l)
                    sum += lower(k, l) * column[l];

                column[k] = (static_cast<T>(j == k) - sum) / lower(k, k);
            }

            lower_inv.column(j, column);
        }

        // Solution.    
        Vector<T> solution{size, 1.0}, old_solution{size};

        // Gauss-Seidel.
        M T_matrix = - (lower_inv * upper);
        Vector<T> C_vector = lower_inv * vector;

        // Residual.
        Vector<T> residual = C_vector;

        // Method.
        do {

            // Solution evaluation.
            old_solution = solution;
            solution = T_matrix * solution + C_vector;

            // Residual evaluation.
            residual = solution - old_solution;

            ++iterations;
        } while((norm(residual) > ALGEBRA_TOLERANCE) || (iterations > ALGEBRA_ITER_MAX));

        #ifdef VERBOSE
        std::cout << "\tConvergence: " << ((iterations >= ALGEBRA_ITER_MAX) ? "failure" : "success") << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << norm(residual) << std::endl;
        #endif

        return solution;
    }
}

#endif