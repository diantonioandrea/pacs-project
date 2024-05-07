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
     * @param a 
     * @param b 
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
        Vector<double> evaluation{compute};

        for(std::size_t j = 0; j < compute; ++j)
            evaluation[j] = std::cos(M_PI * (static_cast<double>(j) + 0.75) / (static_cast<double>(order) + 0.5));

        // std::cout << evaluation << std::endl;

        // Error.
        Vector<double> error_vector{compute, 1.0 + QUADRATURE_TOLERANCE};
        double error = 1.0 + QUADRATURE_TOLERANCE;

        // Utilities.
        Matrix<double> matrix{compute, 3};
        Vector<double> vector{compute};

        // Algorithm.
        while(error > QUADRATURE_TOLERANCE) {

            // Iteration matrix.
            matrix.column(0, 1.0);

            // Iteration.
            for(std::size_t j = 1; j <= order; ++j) {
                double first = static_cast<double>(2 * j - 1) / static_cast<double>(j);
                double second = static_cast<double>(j - 1) / static_cast<double>(j);

                matrix.column(2, matrix.column(1));
                matrix.column(1, matrix.column(0));

                matrix.column(0, first * evaluation * matrix.column(1) - second * matrix.column(2));
            }

            // Step evaluation.
            vector = order * (evaluation * matrix.column(0) - matrix.column(1) / (evaluation * evaluation - 1.0));
            Vector<double> old_evaluation{evaluation};

            for(std::size_t j = 0; j < compute; ++j) {
                double coefficient = (error_vector[j] > QUADRATURE_TOLERANCE) ? matrix(j, 0) / vector[j] : 0.0;
                evaluation[j] = old_evaluation[j] - coefficient;
            }

            // Error evaluation.
            error_vector = std::abs(evaluation - old_evaluation);
            error = max(error_vector);
        }

        // Enlarges evaluation and vector.
        Vector<double> full_evaluation = stack(evaluation(0, compute - 1), flip(evaluation));
        Vector<double> full_evaluation_minus = stack(evaluation(0, compute - 1), -flip(evaluation));
        Vector<double> full_vector = stack(vector(0, compute - 1), flip(vector));

        // Evaluates nodes and weights.
        nodes = ((b + a) / 2) - ((b - a) / 2) * full_evaluation_minus;
        weights = (b - a) / ((1 - (full_evaluation * full_evaluation)) * (full_vector * full_vector));

        return std::vector<Vector<double>>{nodes, weights};
    }

}