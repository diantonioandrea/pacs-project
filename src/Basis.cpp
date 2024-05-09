/**
 * @file Basis.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Basis.hpp>

// Legendre.
#include <Legendre.hpp>

// Assertions.
#include <cassert>

namespace pacs {

    /**
     * @brief Returns the nodal evaluation of a given element's basis.
     * 
     * @param mesh 
     * @param index 
     * @param nodes 
     * @return std::array<Vector<Real>, 3> 
     */
    std::array<Matrix<Real>, 3> basis_2d(const Mesh &mesh, const std::size_t &index, const std::array<Vector<Real>, 2> &nodes) {
        // Element information.
        Element element = mesh.elements[index];
        Polygon polygon = mesh.element(index);
        std::size_t degree = element.degree;

        #ifndef NDEBUG // Integrity checks.
        assert(polygon.edges().size() > 3);
        assert(nodes[0].length == nodes[1].length);
        #endif

        // Element's box.
        auto [xy_min, xy_max] = polygon.box();

        Real x_min = xy_min[0], y_min = xy_min[1];
        Real x_max = xy_max[0], y_max = xy_max[1];

        // Jacobian.
        Matrix<Real> jacobian{2, 2};

        jacobian(0, 0) = 0.5 * (x_max - x_min);
        jacobian(1, 1) = 0.5 * (y_max - y_min);

        // Translation.
        Vector<Real> translation{2};

        translation[0] = 0.5 * (x_max + x_min);
        translation[1] = 0.5 * (y_max + y_min);

        // Jacobian's determinant.
        Real jacobian_det = jacobian(0, 0) * jacobian(1, 1);

        // Inverses.
        Matrix<Real> jacobian_inv{2, 2};
        Vector<Real> translation_inv{2};

        jacobian_inv(0, 0) = jacobian(1, 1) / jacobian_det;
        jacobian_inv(0, 1) = -jacobian(0, 1) / jacobian_det;
        jacobian_inv(1, 0) = -jacobian(1, 0) / jacobian_det;
        jacobian_inv(1, 1) = jacobian(0, 0) / jacobian_det;

        translation_inv[0] = (-jacobian(1, 1) * translation[0] + jacobian(0, 1) * translation[1]) / jacobian_det;
        translation_inv[1] = (jacobian(1, 0) * translation[0] - jacobian(0, 0) * translation[1]) / jacobian_det;

        // Reference nodes.
        Vector<Real> nodes_x{nodes[0].length}, nodes_y{nodes[1].length};

        for(std::size_t j = 0; j < nodes[0].length; ++j) {
            Vector<Real> nodes_j{2};

            nodes_j[0] = nodes[0][j];
            nodes_j[1] = nodes[1][j];

            Vector<Real> product = jacobian_inv * nodes_j + translation_inv;

            nodes_x[j] = product[0];
            nodes_y[j] = product[1];
        }
        
        // Polynomials' degrees.
        std::vector<std::size_t> x_degree, y_degree;

        for(std::size_t j = 0; j < degree + 1; ++j)
            for(std::size_t k = 0; k < degree + 1 - j; ++k) {
                x_degree.emplace_back(j);
                y_degree.emplace_back(k);
            }

        // Evaluations.
        Matrix<Real> phi_evaluations{nodes_x.length, x_degree.size()};
        Matrix<Real> gradx_phi_evaluations{nodes_x.length, x_degree.size()};
        Matrix<Real> grady_phi_evaluations{nodes_x.length, x_degree.size()};

        for(std::size_t j = 0; j < x_degree.size(); ++j) {
            Vector<Real> legendre_x = legendre(nodes_x, x_degree[j]);
            Vector<Real> legendre_y = legendre(nodes_y, y_degree[j]);

            Vector<Real> grad_legendre_x = grad_legendre(nodes_x, x_degree[j]);
            Vector<Real> grad_legendre_y = grad_legendre(nodes_y, y_degree[j]);

            Real coefficient = std::sqrt((2.0 * x_degree[j] + 1.0) * (2.0 * y_degree[j] + 1.0)) / 2.0;

            Vector<Real> phi = coefficient * legendre_x * legendre_y;
            Vector<Real> gradx_phi = coefficient * grad_legendre_x * legendre_y;
            Vector<Real> grady_phi = coefficient * legendre_x * grad_legendre_y;

            phi_evaluations.column(j, phi);
            gradx_phi_evaluations.column(j, gradx_phi);
            grady_phi_evaluations.column(j, grady_phi);
        }

        return {phi_evaluations, gradx_phi_evaluations, grady_phi_evaluations};
    }

}