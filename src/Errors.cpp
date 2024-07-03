/**
 * @file Errors.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Algebra.hpp>
#include <Fem.hpp>
#include <Laplacian.hpp>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Error structure.
     * 
     * @param mesh 
     * @param matrices 
     * @param numerical 
     * @param exact 
     */
    Error::Error(const Mesh &mesh, const std::array<Sparse<Real>, 2> &matrices, const Vector<Real> &numerical, const Functor &exact):
    elements{mesh.elements.size()}, l2_errors{mesh.elements.size()}, h1_errors{mesh.elements.size()} {

        #ifndef NVERBOSE
        std::cout << "Evaluating errors." << std::endl;
        #endif

        // Matrices.
        auto [mass, dg_laplacian] = matrices;

        // Mass blocks.
        auto blocks = block_mass(mesh);

        // Error vector.
        Vector<Real> u_modals = modal(mesh, exact);
        Vector<Real> u_coeff = solve(mass, u_modals, blocks, DB);

        Vector<Real> error = u_coeff - numerical;

        // DG Error.
        this->dg_error = std::sqrt(dot(error, dg_laplacian * error));

        // L2 Error.
        this->l2_error = std::sqrt(dot(error, mass * error));

        // Dofs.
        this->dofs = mesh.dofs();

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(mesh.quadrature);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(mesh.quadrature);

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
                auto [phi, gradx_phi, grady_phi] = basis_2d(mesh, j, {physical_x, physical_y});

                // Solutions.
                Vector<Real> u = exact(physical_x, physical_y);
                Vector<Real> uh = phi * numerical(indices);

                Vector<Real> grad_u = (gradx_phi + grady_phi) * u_coeff(indices);
                Vector<Real> grad_uh = (gradx_phi + grady_phi) * numerical(indices);

                // Local L2 error.
                this->l2_errors[j] += dot(scaled, (u - uh) * (u - uh));

                // Local H1 error.
                this->h1_errors[j] += dot(scaled, (grad_u - grad_uh) * (grad_u - grad_uh));
            }

            this->l2_errors[j] = std::sqrt(this->l2_errors[j]);
            this->h1_errors[j] = std::sqrt(this->h1_errors[j]);

        }

        // Data.
        this->degree = 0;

        for(const auto &element: mesh.elements)
            this->degree = (element.degree > this->degree) ? element.degree : this->degree;

        this->size = 0.0;
        this->elements = mesh.elements.size();

        for(const auto &element: mesh.elements)
            for(const auto &p: element.element.points)
                for(const auto &q: element.element.points)
                    this->size = (distance(p, q) > this->size) ? distance(p, q) : this->size;
    }

    // OUTPUT.
    
    /**
     * @brief Error output.
     * 
     * @param ost 
     * @param error 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Error &error) {
        ost << "Elements: " << error.elements << "\n";
        ost << "Dofs: " << error.dofs << "\n";
        ost << "Degree (p): " << error.degree << "\n";
        ost << "Size (h): " << error.size << "\n";
        ost << "L2 Error: " << error.l2_error << "\n";
        return ost << "DG Error: " << error.dg_error << std::endl;
    }

}