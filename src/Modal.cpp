/**
 * @file Modal.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>
#include <Algebra.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace pacs {

    /**
     * @brief Returns the modal coefficients of a function.
     * 
     * @param mesh 
     * @param function 
     * @return Vector<Real> 
     */
    Vector<Real> modal(const Mesh &mesh, const Functor &function) {

        // Number of quadrature nodes.
        std::size_t degree = mesh.quadrature;

        // Coefficients.
        Vector<Real> coefficients{mesh.dofs()};

        // Quadrature nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {

            // Local dofs.
            std::size_t element_dofs = mesh.elements[j].dofs();

            // Global matrix indices.
            std::vector<std::size_t> indices;

            for(std::size_t k = 0; k < element_dofs; ++k)
                indices.emplace_back(starts[j] + k);
            
            // Polygon.
            Polygon polygon = mesh.element(j);

            // Element sub-triangulation.
            std::vector<Polygon> triangles = triangulate(polygon);

            // Local coefficients.
            Vector<Real> local_coefficients{element_dofs};

            // Loop over the sub-triangulation.
            for(std::size_t k = 0; k < triangles.size(); ++k) {

                // Triangle.
                Polygon triangle = triangles[k];

                // Jacobian.
                Matrix<Real> jacobian{2, 2};

                jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
                jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
                jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
                jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

                // Jacobian's determinant.
                Real jacobian_det = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

                // Translation.
                Vector<Real> translation{2};

                translation[0] = triangle.points[0][0];
                translation[1] = triangle.points[0][1];

                // Physical nodes.
                Vector<Real> physical_x{nodes_x_2d.length};
                Vector<Real> physical_y{nodes_y_2d.length};

                for(std::size_t l = 0; l < physical_x.length; ++l) {
                    Vector<Real> node{2};

                    node[0] = nodes_x_2d[l];
                    node[1] = nodes_y_2d[l];

                    Vector<Real> transformed = jacobian * node + translation;

                    physical_x[l] = transformed[0];
                    physical_y[l] = transformed[1];
                }

                // Weights scaling.
                Vector<Real> scaled = jacobian_det * weights_2d;

                // Basis functions.
                Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
                Matrix<Real> scaled_phi = phi;

                for(std::size_t l = 0; l < scaled_phi.columns; ++l)
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);

                // function solution.
                Vector<Real> local_function = function(physical_x, physical_y);

                // Local coefficients.
                local_coefficients += scaled_phi.transpose() * local_function;
            }

            // Update.
            coefficients(indices, local_coefficients);
        }

        return coefficients;
    }

}