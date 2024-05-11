/**
 * @file Solution.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Solution.hpp>

// Quadrature.
#include <Quadrature.hpp>

// Matrices.
#include <Matrix.hpp>

// Basis.
#include <Basis.hpp>

namespace pacs {

    /**
     * @brief Returns the quadrature nodes evaluation of a given numerical solution and of the exact solution.
     * 
     * @param mesh 
     * @param values 
     * @param exact 
     * @return Vector<Real> 
     */
    std::array<Vector<Real>, 4> solution(const Mesh &mesh, const Vector<Real> &values, const Functor &exact) {

        // Number of quadrature nodes.
        std::size_t degree = (mesh.degree() % 2) ? mesh.degree() : mesh.degree() + 1;

        // Edges.
        std::size_t edges = 0;

        for(std::size_t j = 0; j < mesh.elements_number(); ++j)
            edges += mesh.elements[j].edges.size();

        // Entries.
        std::size_t entries = edges * degree * degree;

        // Return values.
        Vector<Real> x{entries};
        Vector<Real> y{entries};
        Vector<Real> numerical_solution{entries};
        Vector<Real> exact_solution{entries};

        // Quadrature nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Local vectors indices.
        std::vector<std::size_t> local_indices;

        for(std::size_t h = 0; h < degree * degree; ++h)
            local_indices.emplace_back(h);

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

                // Basis functions.
                Matrix<Real> phi = basis_2d(mesh, j, {physical_x, physical_y})[0];

                // Numerical solution.
                Vector<Real> local_numerical = phi * values(indices);

                // Exact solution.
                Vector<Real> local_exact = exact(physical_x, physical_y);

                // Writing.
                x(local_indices, physical_x);
                y(local_indices, physical_y);
                numerical_solution(local_indices, local_numerical);
                exact_solution(local_indices, local_exact);

                // Local indices update.
                for(auto &index: local_indices)
                    index += degree * degree;
            }
        }

        return {x, y, numerical_solution, exact_solution};
    }

}