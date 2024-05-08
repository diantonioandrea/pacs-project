/**
 * @file Laplacian.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Laplacian.hpp>

// Quadrature.
#include <Quadrature.hpp>

// Matrices.
#include <Matrix.hpp>

// Basis functions.
#include <Basis.hpp>

namespace pacs {

    /**
     * @brief Returns the matrix for the Laplacian operator.
     * 
     * @param mesh 
     * @return Sparse<double> 
     */
    Sparse<double> laplacian(const Mesh &mesh) {

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(5);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(5);

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Matrices.
        Sparse<double> A{dofs, dofs};
        
        // To be added.
        // Sparse<double> IA{dofs, dofs};
        // Sparse<double> SA{dofs, dofs};

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements_number(); ++j) {

            // Local dofs.
            std::size_t element_dofs = mesh.elements[j].dofs();

            // Global matrix indices.
            std::vector<std::size_t> indices;

            for(std::size_t k = 0; k < element_dofs; ++k)
                indices.emplace_back(j * element_dofs + k);

            // Polygon.
            Polygon polygon = mesh.element(j);

            // Element sub-triangulation.
            std::vector<Polygon> triangles = triangulate(polygon);

            // Local A.
            Matrix<double> local_A{element_dofs, element_dofs};

            // Loop over the sub-triangulation.
            for(std::size_t k = 0; k < triangles.size(); ++k) {

                // Triangle.
                Polygon triangle = triangles[k];

                // Jacobian.
                Matrix<double> jacobian{2, 2};

                jacobian(0, 0) = triangle.points[1][0] - triangle.points[0][0];
                jacobian(0, 1) = triangle.points[2][0] - triangle.points[0][0];
                jacobian(1, 0) = triangle.points[1][1] - triangle.points[0][1];
                jacobian(1, 1) = triangle.points[2][1] - triangle.points[0][1];

                // Jacobian's determinant.
                double jacobian_det = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

                // Translation.
                Vector<double> translation{2};

                translation[0] = triangle.points[0][0];
                translation[1] = triangle.points[0][1];

                // Physical nodes.
                Vector<double> physical_x{nodes_x_2d.length};
                Vector<double> physical_y{nodes_y_2d.length};

                for(std::size_t l = 0; l < physical_x.length; ++l) {
                    Vector<double> node{2};

                    node[0] = nodes_x_2d[l];
                    node[1] = nodes_y_2d[l];

                    Vector<double> transformed = jacobian * node + translation;

                    physical_x[l] = transformed[0];
                    physical_y[l] = transformed[1];
                }

                // Weights scaling.
                Vector<double> scaled = jacobian_det * weights_2d;

                // Basis function.
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Local matrix assembly.
                Matrix<double> scaled_gradx = gradx_phi;
                Matrix<double> scaled_grady = grady_phi;

                for(std::size_t l = 0; l < scaled_gradx.columns; ++l) {
                    scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
                    scaled_grady.column(l, scaled_grady.column(l) * scaled);
                }

                local_A = local_A + scaled_gradx.transpose() * gradx_phi + scaled_grady.transpose() * grady_phi;
            }

            // Global matrix assembly.
            A.insert(indices, indices, local_A);
        }

        return A;
    }

}