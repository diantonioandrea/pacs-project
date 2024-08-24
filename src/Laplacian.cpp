/**
 * @file Laplacian.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Pacs.hpp>

namespace pacs {

    /**
     * @brief Returns the matrix for the Laplacian operator.
     * 
     * @param mesh 
     * @return Sparse<Real> 
     */
    std::array<Sparse<Real>, 3> laplacian(const Mesh &mesh, const Real &penalty_coefficient) {

        #ifndef NVERBOSE
        std::cout << "Computing the laplacian matrix." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::size_t degree = GAUSS_ORDER;

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

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Volume integrals.

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
                local_M += scaled_phi.transpose() * phi;
                local_A += scaled_gradx.transpose() * gradx_phi + scaled_grady.transpose() * grady_phi;
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
            std::vector<Matrix<Real>> local_IAN;
            std::vector<Matrix<Real>> local_SAN;

            // Penalties.
            Vector<Real> penalties = penalty(mesh, j, penalty_coefficient);

            // Edges.
            std::vector<Segment> edges{polygon.edges()};

            // Loop over faces.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {

                // Neighbour information.
                auto [edge, neighbour, n_edge] = element_neighbours[k];

                // Edge geometry.
                Segment segment{edges[k]};

                // Edge's normal.
                Vector<Real> edge_vector{2};

                edge_vector[0] = segment[1][0] - segment[0][0];
                edge_vector[1] = segment[1][1] - segment[0][1];

                Vector<Real> normal_vector{2};

                normal_vector[0] = edge_vector[1];
                normal_vector[1] = -edge_vector[0];

                normal_vector /= norm(normal_vector);

                // Jacobian.
                Matrix<Real> jacobian{2, 2};

                jacobian(0, 0) = segment[1][0] - segment[0][0];
                jacobian(0, 1) = 0.5 * (segment[1][0] - segment[0][0]);
                jacobian(1, 0) = segment[1][1] - segment[0][1];
                jacobian(1, 1) = 0.5 * (segment[1][1] - segment[0][1]);

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

                Matrix<Real> scaled_grad = normal_vector[0] * scaled_gradx + normal_vector[1] * scaled_grady;

                if(neighbour == -1) { // Boundary edge.

                    local_IA += scaled_grad.transpose() * phi;
                    local_SA += (penalties[k] * scaled_phi).transpose() * phi;

                    // Empty small matrices.
                    local_IAN.emplace_back(Matrix<Real>{1, 1});
                    local_SAN.emplace_back(Matrix<Real>{1, 1});

                } else {

                    local_IA += 0.5 * scaled_grad.transpose() * phi;
                    local_SA += (penalties[k] * scaled_phi).transpose() * phi;

                    // Neighbour's basis function.
                    Matrix<Real> n_phi = basis_2d(mesh, neighbour, {physical_x, physical_y})[0];

                    // Neighbour's local matrix.
                    local_IAN.emplace_back(- 0.5 * scaled_grad.transpose() * n_phi);
                    local_SAN.emplace_back(- (penalties[k] * scaled_phi).transpose() * n_phi);
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
                std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

                for(std::size_t h = 0; h < n_dofs; ++h)
                    n_indices.emplace_back(starts[n_index] + h);

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

    /**
     * @brief Extrapolates blocks (indices) based on mass structure.
     * 
     * @param mesh 
     * @return std::vector<std::array<std::vector<std::size_t>, 2>> 
     */
    std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &mesh) {

        #ifndef NVERBOSE
        std::cout << "Evaluating the mass blocks." << std::endl;
        #endif

        // Elements.
        std::size_t elements = mesh.elements.size();

        // Blocks.
        std::vector<std::array<std::vector<std::size_t>, 2>> blocks;
        blocks.resize(elements);
        std::size_t start = 0;

        // Precomputing dofs.
        std::vector<std::size_t> dofs;
        dofs.resize(elements);

        for(std::size_t j = 0; j < elements; ++j)
            dofs[j] = mesh.elements[j].dofs();

        // Evaluating blocks.
        for(std::size_t j = 0; j < elements; ++j) {
            std::vector<std::size_t> indices;
            indices.resize(dofs[j]);

            for(std::size_t k = 0; k < dofs[j]; ++k)
                indices[k] = start + k;
            
            blocks[j] = std::array<std::vector<std::size_t>, 2>{indices, indices};
            start += dofs[j];
        }

        return blocks;
    }

}