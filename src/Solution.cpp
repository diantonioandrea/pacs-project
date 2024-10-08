/**
 * @file Solution.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace pacs {

    /**
     * @brief Constructs a new Solution structure.
     * 
     * @param mesh Mesh.
     * @param numerical Numerical solution.
     * @param exact Exact solution.
     */
    Solution::Solution(const Mesh &mesh, const Vector<Real> &numerical, const Functor &exact):
    x{mesh.entries}, y{mesh.entries}, numerical{mesh.entries}, exact{mesh.entries} {

        // Number of quadrature nodes.
        std::size_t degree = GAUSS_ORDER;

        // Quadrature nodes.
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Starting indices.
        std::vector<std::size_t> starts;
        starts.emplace_back(0);

        for(std::size_t j = 1; j < mesh.elements.size(); ++j)
            starts.emplace_back(starts[j - 1] + mesh.elements[j - 1].dofs());

        // Local vectors indices.
        std::vector<std::size_t> local_indices;

        for(std::size_t h = 0; h < degree * degree; ++h)
            local_indices.emplace_back(h);

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
                Vector<Real> local_numerical = phi * numerical(indices);

                // Exact solution.
                Vector<Real> local_exact = exact(physical_x, physical_y);

                // Writing.
                this->x(local_indices, physical_x);
                this->y(local_indices, physical_y);
                this->numerical(local_indices, local_numerical);
                this->exact(local_indices, local_exact);

                // Local indices update.
                for(auto &index: local_indices)
                    index += degree * degree;
            }
        }
    }

    // OUTPUT.

    /**
     * @brief Outputs the solution to a file.
     * 
     * @param filename Filename.
     */
    void Solution::write(const std::string &filename) {
        // File loading.
        std::ofstream file{filename};

        file << "@ " << filename << "\n";
        file << "@ contplot.py readable mesh\n";
        file << "@ Structure: [x, y, numerical(x, y), exact(x, y)].\n";

        for(std::size_t j = 0; j < this->x.length; ++j) {
            file << this->x[j] << ",";
            file << this->y[j] << ",";
            file << this->numerical[j] << ",";
            file << this->exact[j] << "\n";
        }
    }

}