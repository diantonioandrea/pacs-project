/**
 * @file Quadrature.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-05-06
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <Algebra.hpp>
#include <Fem.hpp>

#include <cmath>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace pacs {

    /**
     * @brief Returns the Gauss-Legendre nodes and weights of a given order over a given interval [a, b].
     * 
     * @param a 
     * @param b 
     * @param order 
     * @return std::vector<Vector<Real>> 
     */
    std::pair<Vector<Real>, Vector<Real>> gauss_legendre(const Real &a, const Real &b, const std::size_t &order) {
        #ifndef NDEBUG // Integrity check.
        assert(b > a);
        assert(order % 2);
        #endif

        // Nodes and weights.
        Vector<Real> nodes{order}, weights{order};

        // Nodes to be computed before reflection.
        std::size_t compute = (order + 1) / 2;

        // Function evaluation.
        Vector<Real> evaluation{compute};

        for(std::size_t j = 0; j < compute; ++j)
            evaluation[j] = std::cos(M_PI * (static_cast<Real>(j) + 0.75) / (static_cast<Real>(order) + 0.5));

        // std::cout << evaluation << std::endl;

        // Error.
        Vector<Real> error_vector{compute, 1.0 + QUADRATURE_TOLERANCE};
        Real error = 1.0 + QUADRATURE_TOLERANCE;

        // Utilities.
        Matrix<Real> matrix{compute, 3};
        Vector<Real> vector{compute};

        // Algorithm.
        while(error > QUADRATURE_TOLERANCE) {

            // Iteration matrix.
            matrix.column(0, 1.0);

            // Iteration.
            for(std::size_t j = 1; j <= order; ++j) {
                Real first = static_cast<Real>(2 * j - 1) / static_cast<Real>(j);
                Real second = static_cast<Real>(j - 1) / static_cast<Real>(j);

                matrix.column(2, matrix.column(1));
                matrix.column(1, matrix.column(0));

                matrix.column(0, first * evaluation * matrix.column(1) - second * matrix.column(2));
            }

            // Step evaluation.
            vector = order * (evaluation * matrix.column(0) - matrix.column(1) / (evaluation * evaluation - 1.0));
            Vector<Real> old_evaluation{evaluation};

            for(std::size_t j = 0; j < compute; ++j) {
                Real coefficient = (error_vector[j] > QUADRATURE_TOLERANCE) ? matrix(j, 0) / vector[j] : 0.0;
                evaluation[j] = old_evaluation[j] - coefficient;
            }

            // Error evaluation.
            error_vector = std::abs(evaluation - old_evaluation);
            error = max(error_vector);
        }

        // Enlarges evaluation and vector.
        Vector<Real> full_evaluation = stack(evaluation(0, compute - 1), flip(evaluation));
        Vector<Real> full_evaluation_minus = stack(evaluation(0, compute - 1), -flip(evaluation));
        Vector<Real> full_vector = stack(vector(0, compute - 1), flip(vector));

        // Evaluates nodes and weights.
        nodes = ((b + a) / 2) - ((b - a) / 2) * full_evaluation_minus;
        weights = (b - a) / ((1 - (full_evaluation * full_evaluation)) * (full_vector * full_vector));

        return {nodes, weights};
    }

    /**
     * @brief Returns the Gauss-Legendre quadrature nodes and weights of a given order over the reference interval [0, 1].
     * 
     * @param order 
     * @return std::vector<Vector<Real>> 
     */
    std::pair<Vector<Real>, Vector<Real>> quadrature_1d(const std::size_t &order) {

        // Gauss-Legendre nodes over [0, 1].
        return gauss_legendre(0.0, 1.0, order);
    }

    /**
     * @brief Returns the Gauss-Legendre quadrature nodes and weights of a given order over the reference set [-1, 1] x [-1, 1].
     * 
     * @param order 
     * @return std::array<Vector<Real>, 3> 
     */
    std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &order) {

        // Gauss-Legendre nodes over [-1, 1].
        auto [partial_nodes, partial_weights] = gauss_legendre(-1.0, 1.0, order);

        Matrix<Real> row_nodes{order, order}, column_nodes{order, order}, row_weights{order, order}, column_weights{order, order};

        for(std::size_t j = 0; j < order; ++j) {
            row_nodes.row(j, partial_nodes);
            column_nodes.column(j, partial_nodes);
            row_weights.row(j, partial_weights);
            column_weights.column(j, partial_weights);
        }

        Vector<Real> squashed_row_nodes{squash(row_nodes)}, squashed_column_nodes{squash(column_nodes)};
        Vector<Real> squashed_row_weights{squash(row_weights)}, squashed_column_weights{squash(column_weights)};

        // Nodes.
        Vector<Real> nodes_x = (1 + squashed_column_nodes) / 2;
        Vector<Real> nodes_y = (1 - squashed_column_nodes) * (1 + squashed_row_nodes) / 4;

        // Weights.
        Vector<Real> weights = (1 - squashed_column_nodes) * squashed_column_weights * squashed_row_weights / 8;

        return {nodes_x, nodes_y, weights};
    }
    
}