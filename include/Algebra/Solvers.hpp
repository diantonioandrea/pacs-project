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

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 25E3
#endif

// Algebra m limit.
#ifndef ALGEBRA_M_MAX
#define ALGEBRA_M_MAX 25E1
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
     * CGM: Conjugate Gradient method.
     * BICGSTAB: Biconjugate Gradient Stabilized method.
     * 
     */
    enum SparseSolver {GMRES, CGM, BICGSTAB};

    /**
     * @brief Sparse direct solvers.
     * QRB: QR block decomposition method.
     * DB: Diagonal Block method.
     * 
     */
    enum SparseDSolver {DB};

    /**
     * @brief Block sparse preconditioners.
     * DBI: Diagonal Block Inverse method.
     * 
     */
    enum Preconditioner {DBI};

    /**
     * @brief Directly solves a linear system Ax = b.
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
     * @brief Directly solves a linear system Ax = B.
     * 
     * @tparam T 
     * @param A 
     * @param B 
     * @return Matrix<T> 
     */
    template<NumericType T>
    Matrix<T> solve(const Matrix<T> &A, const Matrix<T> &B, const DenseSolver &S = QRD) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == B.rows);
        #endif

        if(S == QRD)
            return _qr(A, B);

        // Default.
        return _qr(A, B);
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
    Vector<T> solve(const Sparse<T> &A, const Vector<T> &b, const SparseSolver &S = GMRES, const Real &TOL = 1E-8) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == b.length);
        #endif

        if(S == GMRES)
            return _gmres(A, b, TOL);

        if(S == CGM)
            return _cgm(A, b, TOL);

        if(S == BICGSTAB)
            return _bicgstab(A, b, TOL);

        // Default.
        return _bicgstab(A, b, TOL);
    }

    /**
     * @brief Solves a sparse linear system MAx = Mb.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param blocks 
     * @param S 
     * @param P 
     * @param TOL 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> solve(const Sparse<T> &A, const Vector<T> &b, std::vector<std::array<std::vector<std::size_t>, 2>> &blocks, const SparseSolver &S = GMRES, const Preconditioner &P = DBI, const Real &TOL = 1E-8) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == b.length);
        #endif

        if(P == DBI) {
            Sparse<T> M = _dbi(A, blocks);
            return solve(M * A, M * b, S, TOL);
        }

        // Default, no preconditioner.
        return solve(A, b, S, TOL);
    }

    /**
     * @brief Directly solves a sparse linear system Ax = b.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param blocks 
     * @param S 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> solve(const Sparse<T> &A, const Vector<T> &b, std::vector<std::array<std::vector<std::size_t>, 2>> &blocks, const SparseDSolver &S = DB) {
        #ifndef NDEBUG // Integrity check.
        assert(A.rows == b.length);
        #endif

        if(S == DB)
            return _db(A, b, blocks);

        // Default.
        return _db(A, b, blocks);
    }

    // MATRIX SOLVERS.

    /**
     * @brief LU solver. No verbosity.
     * 
     * @param vector 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _lu(const Matrix<T> &A, const Vector<T> &b) {

        // LU decomposition.
        auto [L, U] = LU(A);

        // Solves Ly = b using forward substitution.
        Vector<T> y{A.columns};

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
                sum += U(i - 1, j) * x.elements[j];

            x.elements[i - 1] = (y[i - 1] - sum) / U(i - 1, i - 1);
        }

        return x;
    }
    
    /**
     * @brief QR solver. No verbosity.
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

        for (std::size_t i = R.columns; i > 0; --i) {
            T sum = static_cast<T>(0);

            for (std::size_t j = i; j < R.columns; ++j)
                sum += R(i - 1, j) * x.elements[j];

            x.elements[i - 1] = (Qb[i - 1] - sum) / R(i - 1, i - 1);
        }

        return x;
    }

    /**
     * @brief QR Matrix solver. No verbosity.
     * 
     * @tparam T 
     * @param A 
     * @param B 
     * @return Matrix<T> 
     */
    template<NumericType T>
    Matrix<T> _qr(const Matrix<T> &A, const Matrix<T> &B) {
        
        // QR decomposition.
        auto [Q, R] = QR(A);

        // Solution.
        Matrix<T> X{A.columns, B.columns};

        // Solves Rx = QTb for every b using backward substitution.
        for(std::size_t k = 0; k < B.columns; ++k) {
            Vector<T> x{A.columns};
            Vector<T> Qb{Q.transpose() * B.column(k)};

            for (std::size_t i = R.columns; i > 0; --i) {
                T sum = static_cast<T>(0);

                for (std::size_t j = i; j < R.columns; ++j)
                    sum += R(i - 1, j) * x.elements[j];

                x.elements[i - 1] = (Qb[i - 1] - sum) / R(i - 1, i - 1);
            }

            // Update X.
            X.column(k, x);
        }

        return X;
    }

    // SPARSE SOLVERS.

    /**
     * @brief Guessless restarted GMRES.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param TOL 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _gmres(const Sparse<T> &A, const Vector<T> &b, const Real &TOL = 1E-8) {
        return _gmres(A, b, Vector<T>{A.rows}, TOL);
    }

    /**
     * @brief Restarted GMRES.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param guess 
     * @param TOL 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _gmres(const Sparse<T> &A, const Vector<T> &b, const Vector<T> &guess, const Real &TOL = 1E-8) {
        #ifndef NDEBUG
        assert(A.rows == A.columns);
        assert(A.rows == b.length);
        #endif

        #ifndef NVERBOSE
        std::cout << "Solving a linear system with GMRES." << std::endl;
        #endif

        // Problem's size.
        const std::size_t size = A.rows;

        // Iterations.
        std::size_t iterations = 0;

        // m.
        std::size_t m = 1;

        // Solution.
        Vector<T> x{guess};

        // Residual.
        Vector<T> residual = b - A * x;

        #ifndef NVERBOSE
        std::cout << "\tStarting residual: " << norm(residual) << std::endl;
        #endif

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
                    H.elements[k * m + j] = dot(w, Vs[k]);
                    w -= H.elements[k * m + j] * Vs[k];
                }

                // New element for H.
                H.elements[(j + 1) * m + j] = norm(w);

                // New v.
                Vs.emplace_back(w / H.elements[(j + 1) * m + j]);
            }

            // V.
            Matrix<T> V{size, m};

            for(std::size_t j = 0; j < m; ++j)
                for(std::size_t k = 0; k < size; ++k)
                    V.elements[k * m + j] = Vs[j].elements[k];

            // Least squares' right-hand side.
            Vector<T> rhs{m + 1};
            rhs.elements[0] = beta;

            // H rotations.
            Matrix<T> rotation = identity<T>(m + 1);

            for(std::size_t j = 0; j < m; ++j) {
                Matrix<T> rotation_j = rotation;

                // Rotation coefficients.
                T first = H.elements[j * m + j];
                T second = H.elements[(j + 1) * m + j];
                T denominator = std::sqrt(first * first + second * second);

                T s = second / denominator;
                T c = first / denominator;

                rotation_j.elements[j * (m + 1) + j] = c;
                rotation_j.elements[(j + 1) * (m + 1) + j + 1] = c;

                rotation_j.elements[j * (m + 1) + j + 1] = s;
                rotation_j.elements[(j + 1) * (m + 1) + j] = -s;

                // Rotation.
                H = rotation_j * H;
                rhs = rotation_j * rhs;
            }

            // Solves Hy = rhs by backward substitution.
            Vector<T> y{m};

            for(std::size_t j = m; j > 0; --j) {
                #ifdef PARALLEL
                T sum = std::transform_reduce(POLICY, y.elements.begin() + j, y.elements.end(), H.elements.begin() + ((j - 1) * m + j), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * second; });
                #else
                T sum = std::transform_reduce(y.elements.begin() + j, y.elements.end(), H.elements.begin() + ((j - 1) * m + j), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * second; });
                #endif
                
                y.elements[j - 1] = (rhs.elements[j - 1] - sum) / H.elements[(j - 1) * m + j - 1];
            }

            // Solution estimate. 
            x += V * y;

            // Residual update.
            residual = b - A * x;

            // Exit condition.
            if(std::abs(rhs.elements[m]) < TOL)
                break;

            // m update.
            if(m > ALGEBRA_M_MAX) {
                #ifndef NVERBOSE
                std::cout << "\tRestarting, residual: " << norm(b - A * x) << std::endl;
                #endif
                m = 1;
            } else
                ++m;

        } while(iterations < ALGEBRA_ITER_MAX);

        #ifndef NVERBOSE
        std::cout << "Results:" << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << norm(b - A * x) << std::endl;
        #endif

        return x;
    }

    /**
     * @brief Conjugate Gradient.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param TOL 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _cgm(const Sparse<T> &A, const Vector<T> &b, const Real &TOL = 1E-8) {
        #ifndef NDEBUG
        assert(A.rows == A.columns);
        assert(A.rows == b.length);
        #endif

        #ifndef NVERBOSE
        std::cout << "Solving a linear system with CG." << std::endl;
        #endif

        // Iterations.
        std::size_t iterations = 0;

        // Solution.
        Vector<T> x{A.columns};

        // Residual.
        Vector<T> old_residual = b;
        Vector<T> residual = b;

        // Direction.
        Vector<T> old_direction = residual;
        Vector<T> direction = residual;

        // Alpha and Beta.
        T alpha, beta;

        do {
            ++iterations;

            #ifndef NVERBOSE
            if(!(iterations % 500))
                std::cout << "\tCG, iteration: " << iterations << std::endl;
            #endif

            alpha = dot(residual, residual) / dot(direction, A * direction);

            // Solution.
            x += alpha * direction;

            // Residual.
            residual -= alpha * A * direction;

            // Exit condition.
            if(norm(residual) < TOL)
                break;

            // Beta.
            beta = dot(residual, residual) / dot(old_residual, old_residual);

            // Direction.
            direction = residual + beta * direction;

        } while(iterations < ALGEBRA_ITER_MAX);

        #ifndef NVERBOSE
        std::cout << "Results:" << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << norm(b - A * x) << std::endl;
        #endif

        // Fixes missing convergence.
        if(norm(b - A * x) >= TOL) {
            #ifndef NVERBOSE
            std::cout << "Refining." << std::endl;
            #endif
            return _gmres(A, b, x, TOL);
        }

        return x;
    }

    /**
     * @brief Biconjugate Gradient Stabilized method.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param TOL 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _bicgstab(const Sparse<T> &A, const Vector<T> &b, const Real &TOL = 1E-8) {
        #ifndef NDEBUG
        assert(A.rows == A.columns);
        assert(A.rows == b.length);
        #endif

        #ifndef NVERBOSE
        std::cout << "Solving a linear system with BICGSTAB." << std::endl;
        #endif

        // Iterations.
        std::size_t iterations = 0;

        // Solution.
        Vector<T> x{A.columns};

        // Residual.
        Vector<T> old_residual = b;
        Vector<T> residual = b;

        Vector<T> residual_hat = b;

        // Direction.
        Vector<T> direction = residual;

        // Alpha, Beta and Rho.
        T alpha, beta, old_rho = dot(residual_hat, residual), rho;

        do {
            ++iterations;

            #ifndef NVERBOSE
            if(!(iterations % 500))
                std::cout << "\tBICGSTAB, iteration: " << iterations << std::endl;
            #endif

            // Nu.
            Vector<T> nu = A * direction;

            // Alpha.
            alpha = old_rho / dot(residual_hat, nu);

            // h.
            Vector<T> h = x + alpha * direction;
        
            // s.
            Vector<T> s = residual - alpha * nu;

            // Exit condition.
            if(norm(s) < TOL) {
                x = h;
                break;
            }

            // t.
            Vector<T> t = A * s;

            // Omega.
            T omega = dot(t, s) / dot(t, t);

            // x.
            x = h + omega * s;

            // Residual.
            residual = s - omega * t;

            // Exit condition.
            if(norm(residual) < TOL)
                break;

            // Rho.
            rho = dot(residual_hat, residual);

            // Beta.
            beta = (rho / old_rho) * (alpha / omega);
            old_rho = rho;

            // Direction.
            direction = residual + beta * (direction - omega * nu);

        } while(iterations < ALGEBRA_ITER_MAX);

        #ifndef NVERBOSE
        std::cout << "Results:" << std::endl;
        std::cout << "\tIterations: " << iterations << std::endl;
        std::cout << "\tResidual: " << norm(b - A * x) << std::endl;
        #endif

        // Fixes missing convergence.
        if(norm(b - A * x) >= TOL) {
            #ifndef NVERBOSE
            std::cout << "Refining." << std::endl;
            #endif
            return _gmres(A, b, x, TOL);
        }

        return x;
    }

    // SPARSE DIRECT SOLVERS.

    /**
     * @brief Diagonal block method.
     * 
     * @tparam T 
     * @param A 
     * @param b 
     * @param blocks 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> _db(const Sparse<T> &A, const Vector<T> &b, const std::vector<std::array<std::vector<std::size_t>, 2>> &blocks) {

        // Solution.
        Vector<T> x{A.columns};

        #ifndef NVERBOSE
        std::cout << "Solving a linear system with DB." << std::endl;
        #endif

        // Local solutions.
        for(const auto &[rows, columns]: blocks) {
            x(rows, solve(A(rows, columns), b(columns)));
        }

        #ifndef NVERBOSE
        std::cout << "Results:" << std::endl;
        std::cout << "\tResidual: " << norm(b - A * x) << std::endl;
        #endif

        // Possible error.
        if(std::isnan(norm(b - A * x)))
            std::exit(-1);

        return x;
    }

    // BLOCK PRECONDITIONERS.

    template<NumericType T>
    Sparse<T> _dbi(const Sparse<T> &A, const std::vector<std::array<std::vector<std::size_t>, 2>> &blocks) {
        
        // Precondiyioning matrix.
        Sparse<T> M{A.rows, A.columns};

        #ifndef NVERBOSE
        std::cout << "Computing the DBI preconditioner." << std::endl;
        #endif

        // Pre-allocating space for inverses.
        std::vector<Matrix<T>> inverses;

        for(const auto &[rows, columns]: blocks)
            inverses.emplace_back(Matrix<T>{rows.size(), columns.size()});

        // Computing blocks.
        #pragma omp parallel for
        for(std::size_t j = 0; j < blocks.size(); ++j) {

            // Coordinates.
            auto &[rows, columns] = blocks[j];

            #ifndef NDEBUG // Integrity check.
            assert(rows.size() == columns.size());
            #endif

            // Inverting blocks.
            inverses[j] = solve(A(rows, columns), identity<T>(rows.size()));
        }

        // Building M.
        for(std::size_t j = 0; j < blocks.size(); ++j) {
            auto &[rows, columns] = blocks[j];
            M.insert(rows, columns, inverses[j]);
        }

        return M;
    }

}

#endif