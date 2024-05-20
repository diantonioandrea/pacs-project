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

// Penalty.
#include <Penalty.hpp>

namespace pacs {

    /**
     * @brief Returns the matrix for the Laplacian operator.
     * 
     * @param mesh 
     * @return Sparse<Real> 
     */
    std::array<Sparse<Real>, 3> laplacian(const Mesh &mesh) {

        #ifndef NVERBOSE
        std::cout << "Computing the laplacian matrix." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::size_t degree = mesh.quadrature;

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(degree);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

        // Matrices.
        Sparse<Real> M{dofs, dofs};
        Sparse<Real> A{dofs, dofs};
        Sparse<Real> IA{dofs, dofs};
        Sparse<Real> SA{dofs, dofs};

        // Volume integrals.

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

            // Local matrices.
            Matrix<Real> local_M{element_dofs, element_dofs};
            Matrix<Real> local_A{element_dofs, element_dofs};

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
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Some products.
                Matrix<Real> scaled_gradx{gradx_phi};
                Matrix<Real> scaled_grady{grady_phi};
                Matrix<Real> scaled_phi{phi};

                for(std::size_t l = 0; l < scaled_gradx.columns; ++l) {
                    scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
                    scaled_grady.column(l, scaled_grady.column(l) * scaled);
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);
                }

                // Local matrix assembly.
                local_M = local_M + scaled_phi.transpose() * phi;
                local_A = local_A + scaled_gradx.transpose() * gradx_phi + scaled_grady.transpose() * grady_phi;
            }

            // Global matrix assembly.
            M.insert(indices, indices, local_M);
            A.insert(indices, indices, local_A);

            // Face integrals.

            // Local matrices.
            Matrix<Real> local_IA{element_dofs, element_dofs};
            Matrix<Real> local_SA{element_dofs, element_dofs};

            // Element's neighbours.
            std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

            // Local matrices for neighbours.
            std::vector<Matrix<Real>> local_IAN(polygon.edges().size(), Matrix<Real>{element_dofs, element_dofs});
            std::vector<Matrix<Real>> local_SAN(polygon.edges().size(), Matrix<Real>{element_dofs, element_dofs});

            // Penalties.
            Vector<Real> penalties = penalty(mesh, j);

            // Loop over faces.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {

                // Neighbour information.
                auto [edge, neighbour, n_edge] = element_neighbours[k];

                // Edge geometry.
                Segment segment{mesh.edge(mesh.elements[j].edges[k])};

                // Edge's normal. Mind the order.
                Vector<Real> unit_vector{2};

                unit_vector[0] = segment[0][1] - segment[1][1];
                unit_vector[1] = segment[1][0] - segment[0][0];

                unit_vector /= unit_vector.norm();

                // Jacobian.
                Matrix<Real> jacobian{2, 2};

                jacobian(0, 0) = segment[1][0] - segment[0][0];
                jacobian(0, 1) = 0.5 * jacobian(0, 0);
                jacobian(1, 0) = segment[1][1] - segment[0][1];
                jacobian(1, 1) = 0.5 * jacobian(1, 0);

                // Translation.
                Vector<Real> translation{2};

                translation[0] = segment[0][0];
                translation[1] = segment[0][1];

                // Physical nodes.
                Vector<Real> physical_x{nodes_1d.length};
                Vector<Real> physical_y{nodes_1d.length};

                for(std::size_t l = 0; l < nodes_1d.length; ++l) {
                    Vector<Real> node{2};

                    node[0] = nodes_1d[l];
                    node[1] = 0.0;

                    Vector<Real> transformed = jacobian * node + translation;

                    physical_x[l] = transformed[0];
                    physical_y[l] = transformed[1];
                }

                // Weights scaling.
                Vector<Real> scaled = std::abs(segment) * weights_1d;

                // Basis functions.
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Local matrix assembly.
                Matrix<Real> scaled_gradx{gradx_phi};
                Matrix<Real> scaled_grady{grady_phi};
                Matrix<Real> scaled_phi{phi};

                for(std::size_t l = 0; l < scaled_gradx.columns; ++l) {
                    scaled_gradx.column(l, scaled_gradx.column(l) * scaled);
                    scaled_grady.column(l, scaled_grady.column(l) * scaled);
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);
                }

                Matrix<Real> scaled_grad = unit_vector[0] * scaled_gradx + unit_vector[1] * scaled_grady;

                if(neighbour == -1) { // Boundary edge.

                    local_IA = local_IA + scaled_grad.transpose() * phi;
                    local_SA = local_SA + (penalties[k] * scaled_phi).transpose() * phi;

                } else {

                    local_IA = local_IA + 0.5 * scaled_grad.transpose() * phi;
                    local_SA = local_SA + (penalties[k] * scaled_phi).transpose() * phi;

                    // Neighbour's basis function.
                    Matrix<Real> n_phi = basis_2d(mesh, neighbour, {physical_x, physical_y})[0];

                    // Neighbour's local matrix.
                    local_IAN[k] = local_IAN[k] - 0.5 * scaled_grad.transpose() * n_phi;
                    local_SAN[k] = local_SAN[k] - (penalties[k] * scaled_phi).transpose() * n_phi;
                }
            }

            IA.insert(indices, indices, local_IA);
            SA.insert(indices, indices, local_SA);

            // Neighbouring DG matrices assembly.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {
                if(element_neighbours[k][1] == -1)
                    continue;

                std::vector<std::size_t> n_indices;
                std::size_t n_index = element_neighbours[k][1];
                std::size_t n_dofs = mesh.elements[n_index].dofs();

                for(std::size_t h = 0; h < n_dofs; ++h)
                    n_indices.emplace_back(n_index * n_dofs + h);

                IA.add(indices, n_indices, local_IAN[k]);
                SA.add(indices, n_indices, local_SAN[k]);
            }
        }

        // Matrices.
        Sparse<Real> mass = M;
        Sparse<Real> dg_stiffness = A + SA;
        Sparse<Real> stiffness = dg_stiffness - IA - IA.transpose();

        // Compression.
        mass.compress();
        dg_stiffness.compress();
        stiffness.compress();
        
        return {mass, stiffness, dg_stiffness};
    }

}