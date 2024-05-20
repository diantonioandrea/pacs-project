/**
 * @file Forcing.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Forcing.hpp>

// Quadrature.
#include <Quadrature.hpp>

// Vectors.
#include <Vector.hpp>

// Matrices.
#include <Matrix.hpp>

// OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// Basis functions.
#include <Basis.hpp>

namespace pacs {

    /**
     * @brief Assemblies the RHS.
     * 
     * @param mesh 
     * @return Vector<Real> 
     */
    Vector<Real> forcing(const Mesh &mesh, const Functor &source) {

        #ifndef NVERBOSE
        std::cout << "Computing the forcing term." << std::endl;
        #endif

        // Number of quadrature nodes.
        std::size_t degree = mesh.quadrature;

        // Quadrature nodes.
        auto [nodes_1d, weights_1d] = quadrature_1d(degree);
        auto [nodes_x_2d, nodes_y_2d, weights_2d] = quadrature_2d(degree);

        // Degrees of freedom.
        std::size_t dofs = mesh.dofs();

        // Forcing term.
        Vector<Real> forcing{dofs};

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

            // Local forcing term.
            Vector<Real> local_f{element_dofs};

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

                // Local source evaluation.
                Vector<Real> local_source = source(physical_x, physical_y);

                // Basis functions.
                auto phi = basis_2d(mesh, j, {physical_x, physical_y})[0];
                Matrix<Real> scaled_phi{phi};

                for(std::size_t l = 0; l < scaled_phi.columns; ++l)
                    scaled_phi.column(l, scaled_phi.column(l) * scaled);

                // Local forcing term.
                local_f = local_f + scaled_phi.transpose() * local_source;
            }

            // Global forcing term.
            forcing(indices, local_f);

        }

        return forcing;

    }

}