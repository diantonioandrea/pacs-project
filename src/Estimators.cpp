/**
 * @file Estimators.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Algebra.hpp>
#include <Fem.hpp>
#include <Laplacian.hpp>

namespace pacs {

    /**
     * @brief Construct a new Estimator structure.
     * 
     * @param mesh 
     * @param matrices 
     * @param numerical 
     * @param forcing 
     * @param source 
     * @param dirichlet 
     */
    Estimator::Estimator(const Mesh &mesh, const std::array<Sparse<Real>, 2> &matrices, const Vector<Real> &numerical, const Functor &source, const Functor &dirichlet):
    elements{mesh.elements.size()}, estimates{mesh.elements.size()} {
        
        #ifndef NVERBOSE
        std::cout << "Evaluating estimates." << std::endl;
        #endif

        // Matrices.
        auto [mass, dg_laplacian] = matrices;

        // Number of quadrature nodes.
        std::size_t degree = mesh.quadrature;

        // Degrees of freedom.
        this->dofs = mesh.dofs();

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(degree);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours = mesh.neighbours;

        // Sizes.
        Vector<Real> sizes{mesh.elements.size()};

        for(std::size_t j = 0; j < sizes.length; ++j) {
            Element element{mesh.elements[j]};

            for(const auto &p: element.element.points)
                for(const auto &q: element.element.points)
                    sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];
        }

        // Loop over the elements.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {

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
                Matrix<Real> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

                // Local numerical laplacian.
                Vector<Real> lap_uh = lap_phi * numerical(indices);

                // Local exact forcing.
                Vector<Real> f = source(physical_x, physical_y);

                // Local estimator, R_{K, E}.
                this->estimates[j] += sizes[j] * sizes[j] * dot(scaled * (f + lap_uh), f + lap_uh);
            }

            // Element's neighbours.
            std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

            // Penalties.
            Vector<Real> penalties = penalty(mesh, j);

            // Edges.
            std::vector<Segment> edges{polygon.edges()};

            // Loop over faces.
            for(std::size_t k = 0; k < element_neighbours.size(); ++k) {

                // Neighbour information.
                auto [edge, neighbour, n_edge] = element_neighbours[k];

                // Edge geometry.
                Segment segment{edges[k]}; // Mesh's edges to be fixed. [!]

                // Edge's normal. Check the order. [!]
                Vector<Real> edge_vector{2};

                edge_vector[0] = segment[1][0] - segment[0][0];
                edge_vector[1] = segment[1][1] - segment[0][1];

                Vector<Real> normal_vector{2};

                normal_vector[0] = edge_vector[1];
                normal_vector[1] = -edge_vector[0];

                normal_vector /= norm(normal_vector);
                edge_vector /= norm(edge_vector);

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

                // Local numerical solution and gradient.
                Vector<Real> uh = phi * numerical(indices);

                Matrix<Real> grad = normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
                Vector<Real> grad_uh = grad * numerical(indices);

                Matrix<Real> grad_t = edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
                Vector<Real> grad_uh_t = grad * numerical(indices);

                if(neighbour == -1) { // Boundary edge.

                    // Local exact Dirichlet.
                    Vector<Real> g = dirichlet(physical_x, physical_y);

                    // Local estimator, R_{K, J}.
                    this->estimates[j] += penalties[k] * dot(scaled * (uh - g), uh - g);

                    // Local estimator, R_{K, T}.
                    // TBA.

                } else {

                    // Neighbour's basis function.
                    auto [n_phi, n_gradx_phi, n_grady_phi] = basis_2d(mesh, neighbour, {physical_x, physical_y});

                    std::vector<std::size_t> n_indices;
                    std::size_t n_index = element_neighbours[k][1];
                    std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

                    for(std::size_t h = 0; h < n_dofs; ++h)
                        n_indices.emplace_back(n_index * n_dofs + h);

                    // Neighbour's numerical solution and gradient.
                    Vector<Real> n_uh = n_phi * numerical(n_indices);

                    Matrix<Real> n_grad = normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
                    Vector<Real> n_grad_uh = grad * numerical(n_indices);

                    Matrix<Real> n_grad_t = edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
                    Vector<Real> n_grad_uh_t = grad * numerical(n_indices);

                    // Local estimator, R_{K, J}.
                    this->estimates[j] += penalties[k] * dot(scaled * (uh - n_uh), uh - n_uh);

                    // Local estimator, R_{K, N}.
                    this->estimates[j] += sizes[j] * dot(scaled * (grad_uh - n_grad_uh), grad_uh - n_grad_uh);

                    // Local estimator, R_{K, T}.
                    this->estimates[j] += sizes[j] * dot(scaled * (grad_uh_t - n_grad_uh_t), grad_uh_t - n_grad_uh_t);
                }
            }

            this->estimate += this->estimates[j];
            this->estimates[j] = std::sqrt(this->estimates[j]);
        }

        this->estimate = std::sqrt(this->estimate);
    }

    /**
     * @brief Estimator output.
     * 
     * @param ost 
     * @param estimator 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Estimator &estimator) {
        ost << "Elements: " << estimator.elements << std::endl;
        ost << "Dofs: " << estimator.dofs << std::endl;
        return ost << "Estimate: " << estimator.estimate << std::endl;
    }

}