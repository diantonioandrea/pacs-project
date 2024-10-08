/**
 * @file Estimators.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

namespace pacs {

    /**
     * @brief Construct a new Estimator structure.
     * 
     * @param mesh Mesh.
     * @param numerical Numerical solution.
     * @param forcing RHS.
     * @param source Source.
     * @param dirichlet Dirichlet BC.
     */
    Estimator::Estimator(const Mesh &mesh, const Sparse<Real> &mass, const Vector<Real> &numerical, const Functor &source, const Functor &dirichlet, const TwoFunctor &dirichlet_gradient, const Real &penalty_coefficient):
    estimates{mesh.elements.size()}, fits{mesh.elements.size()} {
        
        #ifndef NVERBOSE
        std::cout << "Evaluating estimates." << std::endl;
        #endif

        // Degrees of freedom.
        this->dofs = mesh.dofs();

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(GAUSS_ORDER);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(GAUSS_ORDER);

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

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Coefficients.
        Vector<Real> f_modals = modal(mesh, source);
        Vector<Real> g_modals = modal(mesh, dirichlet);

        Vector<Real> f_coeff = (norm(f_modals) > TOLERANCE) ? solve(mass, modal(mesh, source), blocks, DB) : Vector<Real>{mesh.dofs()};
        Vector<Real> g_coeff = (norm(g_modals) > TOLERANCE) ? solve(mass, modal(mesh, dirichlet), blocks, DB) : Vector<Real>{mesh.dofs()};

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
                Matrix<Real> lap_phi = lap_basis_2d(mesh, j, {physical_x, physical_y});

                // Local numerical laplacian.
                Vector<Real> lap_uh = lap_phi * numerical(indices);

                // Local exact source.
                Vector<Real> f = source(physical_x, physical_y);

                // Local source approximation.
                Vector<Real> f_bar = phi * f_coeff(indices);

                // Local estimator, R_{K, E}^2.
                this->estimates[j] += sizes[j] * sizes[j] * dot(scaled, (f_bar + lap_uh) * (f_bar + lap_uh));

                // Local data oscillation, O_{K, E}^2.
                this->estimates[j] += sizes[j] * sizes[j] * dot(scaled, (f - f_bar) * (f - f_bar));
            }

            // Element's neighbours.
            std::vector<std::array<int, 3>> element_neighbours = neighbours[j];

            // Penalties.
            Vector<Real> penalties = penalty(mesh, j, penalty_coefficient);

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
                Matrix<Real> scaled_phi{phi};

                for(std::size_t l = 0; l < scaled_phi.columns; ++l)
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);

                // Local numerical solution and gradients.
                Vector<Real> uh = phi * numerical(indices);

                Matrix<Real> grad = normal_vector[0] * gradx_phi + normal_vector[1] * grady_phi;
                Vector<Real> grad_uh = grad * numerical(indices);

                Matrix<Real> grad_t = edge_vector[0] * gradx_phi + edge_vector[1] * grady_phi;
                Vector<Real> grad_uh_t = grad_t * numerical(indices);

                if(neighbour == -1) { // Boundary edge.

                    // Local exact Dirichlet and gradient.
                    Vector<Real> g = dirichlet(physical_x, physical_y);
                    auto [grad_g_x, grad_g_y] = dirichlet_gradient(physical_x, physical_y);
                    Vector<Real> grad_g_t = edge_vector[0] * grad_g_x + edge_vector[1] * grad_g_y;

                    // Approximate Dirichlet and gradient.
                    Vector<Real> g_bar = phi * g_coeff(indices);
                    Vector<Real> grad_g_t_bar = grad_t * g_coeff(indices);

                    // Local estimator, R_{K, J}^2.
                    this->estimates[j] += penalties[k] * dot(scaled, (uh - g_bar) * (uh - g_bar));

                    // // Local estimator, R_{K, N}^2.
                    // this->estimates[j] += sizes[j] * dot(scaled, grad_uh * grad_uh);

                    // Local estimator, R_{K, T}^2.
                    this->estimates[j] += sizes[j] * dot(scaled, (grad_uh_t - grad_g_t_bar) * (grad_uh_t - grad_g_t_bar));

                    // Local data oscillation, O_{K, J}^2.
                    this->estimates[j] += penalties[k] * dot(scaled, (g - g_bar) * (g - g_bar));

                    // Local data oscillation, O_{K, T}^2.
                    this->estimates[j] += sizes[j] * dot(scaled, (grad_g_t - grad_g_t_bar) * (grad_g_t - grad_g_t_bar));

                } else {

                    // Neighbour's basis function.
                    auto [n_phi, n_gradx_phi, n_grady_phi] = basis_2d(mesh, neighbour, {physical_x, physical_y});

                    std::vector<std::size_t> n_indices;
                    std::size_t n_index = element_neighbours[k][1];
                    std::size_t n_dofs = mesh.elements[n_index].dofs(); // Neighbour's dofs.

                    for(std::size_t h = 0; h < n_dofs; ++h)
                        n_indices.emplace_back(starts[n_index] + h);

                    // Neighbour's numerical solution and gradients.
                    Vector<Real> n_uh = n_phi * numerical(n_indices);

                    Matrix<Real> n_grad = normal_vector[0] * n_gradx_phi + normal_vector[1] * n_grady_phi;
                    Vector<Real> n_grad_uh = n_grad * numerical(n_indices);

                    Matrix<Real> n_grad_t = edge_vector[0] * n_gradx_phi + edge_vector[1] * n_grady_phi;
                    Vector<Real> n_grad_uh_t = n_grad_t * numerical(n_indices);

                    // Local estimator, R_{K, J}^2.
                    this->estimates[j] += penalties[k] * dot(scaled, (uh - n_uh) * (uh - n_uh));

                    // Local estimator, R_{K, N}^2.
                    this->estimates[j] += sizes[j] * dot(scaled, (grad_uh - n_grad_uh) * (grad_uh - n_grad_uh));

                    // Local estimator, R_{K, T}^2.
                    this->estimates[j] += sizes[j] * dot(scaled, (grad_uh_t - n_grad_uh_t) * (grad_uh_t - n_grad_uh_t));
                }
            }

            this->estimate += this->estimates[j];
            this->estimates[j] = std::sqrt(this->estimates[j]);

            // Degrees.
            Vector<Real> degrees{element_dofs};
            std::size_t counter = 0;

            for(std::size_t i = 0; i < mesh.elements[j].degree + 1; ++i)
                for(std::size_t k = 0; k < mesh.elements[j].degree + 1 - i; ++k) {
                    degrees[counter] = static_cast<Real>(i + k);
                    ++counter;
                }
            
            // Coefficients.
            Vector<Real> coefficients = numerical(indices);

            for(auto &coefficient: coefficients.elements)
                coefficient = std::log(std::abs(coefficient));

            // Fit.
            Vector<Real> fit = polyfit(degrees, coefficients, 1);
            this->fits[j] = -fit[1];
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
        ost << "Dofs: " << estimator.dofs << std::endl;
        return ost << "Estimate: " << estimator.estimate << std::endl;
    }

}