/**
 * @file Quadrature.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Quadrature.hpp>

// Math.
#include <cmath>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

// Matrices.
#include <Matrix.hpp>

namespace pacs {

    /**
     * @brief Returns the Gauss-Legendre nodes and weights of a given order over a given interval [a, b].
     * 
     * @param interval 
     * @param order 
     * @return std::vector<Vector<double>> 
     */
    std::vector<Vector<double>> gauss_legendre(const double &a, const double &b, const std::size_t &order) {
        #ifndef NDEBUG // Integrity check.
        assert(b > a);
        assert(order % 2);
        #endif

        // Nodes and weights.
        Vector<double> nodes{order}, weights{order};

        // Nodes to be computed before reflection.
        std::size_t compute = (order + 1) / 2;

        // Function evaluation.
        Vector<double> evaluation{compute}, old_evaluation{compute};

        for(std::size_t j = 0; j < compute; ++j)
            evaluation[j] = std::cos(M_PI * (static_cast<double>(j) + 0.75) / (static_cast<double>(order) + 0.5));

        // Error.
        Vector<double> error_vector{compute, 1.0 + QUADRATURE_TOLERANCE};
        double error = 1.0 + QUADRATURE_TOLERANCE;

        // Utilities.
        Matrix<double> matrix{compute, 3};
        Vector<double> vector{compute};

        // Algorithm.
        while(error > QUADRATURE_TOLERANCE) {
            
            // Iteration matrix.
            for(std::size_t j = 0; j < compute; ++j)
                matrix(j, 0) = 1.0;

            // Iteration.
            for(std::size_t j = 0; j < order; ++j) {
                double first = static_cast<double>(2 * (j + 1) - 1) / static_cast<double>(j + 1);
                double second = static_cast<double>(j) / static_cast<double>(j + 1);

                // Matrix update.
                for(std::size_t k = 0; k < compute; ++k) {
                    matrix(k, 2) = matrix(k, 1);
                    matrix(k, 1) = matrix(k, 0);
                }

                for(std::size_t k = 0; k < compute; ++k) {
                    matrix(k, 0) = first * evaluation[k] * matrix(k, 1) - second * matrix(k, 2);
                }
            }

            // Step evaluation.
            for(std::size_t j = 0; j < compute; ++j) {
                vector[j] = order * (evaluation[j] * matrix(j, 0) - matrix(j, 1)) / (evaluation[j] * evaluation[j] - 1.0);
            }

            old_evaluation = evaluation;

            for(std::size_t j = 0; j < compute; ++j) {
                double coefficient = (error_vector[j] > QUADRATURE_TOLERANCE) ? matrix(j, 0) / vector[j] : 0.0;
                evaluation[j] = old_evaluation[j] - coefficient;
            }

            // Error evaluation.
            error_vector = std::abs(evaluation - old_evaluation);
            error = max(error_vector);
        }

        // Enlarges evaluation and vector.
        Vector<double> long_evaluation{order}, long_flipped_evaluation{order}, long_vector{order};

        for(std::size_t j = 0; j < compute - 1; ++j) {
            long_evaluation[j] = evaluation[j];
            long_flipped_evaluation[j] = evaluation[j];
            long_vector[j] = vector[j];
        }

        for(std::size_t j = 0; j < compute; ++j) {
            long_evaluation[compute + j - 1] = evaluation[compute - j - 1];
            long_flipped_evaluation[compute + j - 1] = -evaluation[compute - j - 1];
            long_vector[compute + j - 1] = vector[compute - j - 1];
        }

        // Nodes computation.
        for(std::size_t j = 0; j < order; ++j)
            nodes[j] = ((b + a) / 2) - ((b - a) / 2) * long_flipped_evaluation[j];

        // Weights evaluation.
        for(std::size_t j = 0; j < order; ++j)
            weights[j] = (b - a) / ((1 - long_evaluation[j] * long_evaluation[j]) * (long_vector[j] * long_vector[j]));

        return {nodes, weights};
    }

}