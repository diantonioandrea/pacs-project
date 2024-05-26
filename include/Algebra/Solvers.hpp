/**
 * @file Solvers.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SOLVERS_PACS
#define SOLVERS_PACS

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Sparse.hpp"

#include "Methods/Vector.hpp"
#include "Methods/Matrix.hpp"
#include "Methods/Sparse.hpp"

// Algebra tolerance.
#ifndef ALGEBRA_TOLERANCE
#define ALGEBRA_TOLERANCE 1E-12
#endif

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 1E3
#endif

// Algebra m limit.
#ifndef ALGEBRA_M_MAX
#define ALGEBRA_M_MAX 2E2
#endif

namespace pacs {

    /**
     * @brief Dense solvers.
     * LUD: LU decomposition method.
     * QRD: QR decomposition method.
     * 
     */
    enum DenseSolver {LUD, QRD};

    /**
     * @brief Sparse solvers.
     * GMRES: Generalized Minimum Residual method.
     * 
     */
    enum SparseSolver {GMRES};

    /**
     * @brief Solves a linear system Ax = b.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param S 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> solve(const Matrix<T> &A, const Vector<T> &b, const DenseSolver &S = QRD) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == b.length);
        #endif

        if(S == LUD)
            return _lu(A, b);

        if(S == QRD)
            return _qr(A, b);

        // Default.
        return _qr(A, b);
    }

    /**
     * @brief Solves a sparse linear system Ax = b.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param S 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> solve(const Sparse<T> &A, const Vector<T> &b, const SparseSolver &S = GMRES) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == b.length);
        #endif

        if(S == GMRES)
            return _gmres(A, b);

        // Default.
        return _gmres(A, b);
    }

    // MATRIX SOLVERS.

    /**
     * @brief LU solver.
     * 
     * @param vector 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _lu(const Matrix<T> &A, const Vector<T> &b) {

        // LU decomposition.
        auto [L, U] = LU(A);

        // Solves Ly = b using forward substitution.
        Vector<T> y{A.rows};

        for (std::size_t i = 0; i < A.rows; ++i) {
            T sum = static_cast<T>(0);

            for (std::size_t j = 0; j < i; ++j)
                sum += L(i, j) * y[j];

            y[i] = (b[i] - sum) / L(i, i);
        }

        // Solves Ux = y using backward substitution.
        Vector<T> x{A.columns};

        for (std::size_t i = A.columns; i > 0; --i) {
            T sum = static_cast<T>(0);

            for (std::size_t j = i; j < A.columns; ++j)
                sum += U(i - 1, j) * x[j];

            x[i - 1] = (y[i - 1] - sum) / U(i - 1, i - 1);
        }

        return x;
    }
    
    /**
     * @brief QR solver.
     * 
     * @param vector 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _qr(const Matrix<T> &A, const Vector<T> &b) {
        
        // QR decomposition.
        auto [Q, R] = QR(A);

        // Solves Rx = QTb using backward substitution.
        Vector<T> x{A.columns};
        Vector<T> Qb{Q.transpose() * b};

        for (std::size_t i = A.columns; i > 0; --i) {
            T sum = static_cast<T>(0);

            for (std::size_t j = i; j < A.columns; ++j)
                sum += R(i - 1, j) * x[j];

            x[i - 1] = (Qb[i - 1] - sum) / R(i - 1, i - 1);
        }

        return x;
    }

    // SPARSE SOLVERS.

    /**
     * @brief Restarted GMRES.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _gmres(const Sparse<T> &A, const Vector<T> &b) {
        #ifndef NDEBUG
        assert(A.rows == A.columns);
        assert(A.rows == b.length);
        #endif

        #ifndef NVERBOSE
        std::cout << "Solving a linear system." << std::endl;
        #endif

        // Problem's size.
        const std::size_t size = A.rows;

        // Iterations.
        std::size_t iterations = 0;

        // m.
        std::size_t m = 1;

        // Solution.
        Vector<T> x{size};

        // Residual.
        Vector<T> residual = b;

        do {
            ++iterations;

            #ifndef NVERBOSE
            if(!(iterations % 50))
                std::cout << "\tRestarted GMRES, iteration: " << iterations << std::endl;
            #endif

            // Beta.
            Real beta = norm(residual);

            // H.
            Matrix<T> H{m + 1, m};

            // Vs.
            std::vector<Vector<T>> Vs;
            Vs.emplace_back(residual / beta);

            // Arnoldi.
            for(std::size_t j = 0; j < m; ++j) {
                Vector<T> w = A * Vs[j];

                for(std::size_t k = 0; k <= j; ++k) {
                    H(k, j) = dot(w, Vs[k]);
                    w -= H(k, j) * Vs[k];
                }

                // New element for H.
                H(j + 1, j) = norm(w);

                // New v.
                Vs.emplace_back(w / H(j + 1, j));
            }

            // V.
            Matrix<T> V{size, m};

            for(std::size_t j = 0; j < m; ++j)
                V.column(j, Vs[j]);

            // Least square's right-hand side.
            Vector<T> rhs{m + 1};
            rhs[0] = beta;

            // H rotations.
            Matrix<T> rotation = identity<T>(m + 1);

            for(std::size_t j = 0; j < m; ++j) {
                Matrix<T> rotation_j = rotation;

                // Rotation coefficients.
                T s = H(j + 1, j) / std::sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));
                T c = H(j, j) / std::sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));

                rotation_j(j, j) = c;
                rotation_j(j + 1, j + 1) = c;

                rotation_j(j, j + 1) = s;
                rotation_j(j + 1, j) = -s;

                // Rotation.
                H = rotation_j * H;
                rhs = rotation_j * rhs;
            }

            // Solves Hy = rhs by backward substitution.
            Vector<T> y{m};

            for(std::size_t j = m; j > 0; --j) {
                T sum = static_cast<T>(0);

                for(std::size_t k = j; k < m; ++k)
                    sum += H(j - 1, k) * y[k];
                
                y[j - 1] = (rhs[j - 1] - sum) / H(j - 1, j - 1);
            }

            // Solution estimate. 
            x = x + V * y;

            // Residual update.
            residual = b - A * x;

            // Exit condition.
            if(std::abs(rhs[m]) < ALGEBRA_TOLERANCE)
                break;

            // m update.
            m = (m <= ALGEBRA_M_MAX) ? m + 1 : 1;

        } while(iterations < ALGEBRA_ITER_MAX);

        #ifndef NVERBOSE
        std::cout << "Results:" << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << norm(b - A * x) << std::endl;
        #endif

        return x;
    }

}

#endif