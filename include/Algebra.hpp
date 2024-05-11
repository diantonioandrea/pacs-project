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
#define ALGEBRA_TOLERANCE 1E-6
#endif

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 2E4
#endif

namespace pacs {

    /**
     * @brief Sparse solver using Gauss-Seidel.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> solve(const Sparse<T> &matrix, const Vector<T> &vector) {
        #ifndef NDEBUG // Integrity checks.
        assert(matrix.rows == matrix.columns);
        assert(matrix.columns == vector.length);
        #endif

        #ifndef NVERBOSE
        std::cout << "Solving a linear system." << std::endl;
        #endif

        // Matrix size.
        const std::size_t size = matrix.rows;

        // Iterations.
        std::size_t iterations = 0;

        // Diagonal's inverse.
        Sparse<T> m_diagonal = matrix.diagonal();
        Sparse<T> diagonal_inv{size, size};

        for(std::size_t j = 0; j < size; ++j)
            diagonal_inv.insert(j, j, static_cast<T>(1) / m_diagonal(j, j));

        // Preconditioned.
        Sparse<T> p_matrix = matrix * diagonal_inv;

        // L, diagonal and U.
        Sparse<T> lower = p_matrix.lower();
        Sparse<T> diagonal = p_matrix.diagonal();
        Sparse<T> upper = p_matrix.upper();

        Sparse<T> LD = lower + diagonal;

        // (L + diagonal)'s determinant.
        T check = mtrace(matrix);
        
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(check) > ALGEBRA_TOLERANCE);
        #endif

        // L's inverse.
        Sparse<T> LD_inv{size, size};

        for(std::size_t j = 0; j < size; ++j) {
            Vector<T> column{size};

            for(std::size_t k = 0; k < size; ++k) {
                T sum = static_cast<T>(0);

                for(std::size_t l = 0; l < k; ++l)
                    sum += LD(k, l) * column[l];

                column[k] = (static_cast<T>(j == k) - sum) / LD(k, k);
            }

            LD_inv.column(j, column);
        }

        // Solution.    
        Vector<T> solution{size, 1.0}, old_solution{size};

        // Gauss-Seidel.
        Sparse<T> T_matrix = - (LD_inv * upper);
        Vector<T> C_vector = LD_inv * vector;

        // Compression.
        T_matrix.compress();

        // Residual.
        Real residual = 1.0;

        // Method.
        do {

            // Solution evaluation.
            old_solution = solution;
            solution = T_matrix * solution + C_vector;

            // Residual evaluation.
            residual = (solution - old_solution).norm();

            ++iterations;
        } while((residual > ALGEBRA_TOLERANCE) || (iterations > ALGEBRA_ITER_MAX));

        #ifndef NVERBOSE
        std::cout << "\tConvergence: " << ((iterations >= ALGEBRA_ITER_MAX) ? "failure" : "success") << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << residual << std::endl;
        #endif

        return diagonal_inv * solution;
    }
}

#endif