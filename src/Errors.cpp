/**
 * @file Errors.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Errors.hpp>

// Modal coefficients.
#include <Solution.hpp>

// Vectors.
#include <Vector.hpp>

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
    l2_errors{mesh.elements.size()} {

        #ifndef NVERBOSE
        std::cout << "Evaluating errors." << std::endl;
        #endif

        // Matrices.
        auto [mass, dg_laplacian] = matrices;

        // Error vector.
        Vector<Real> error = mass.solve<GMRES>(modal(mesh, exact)) - numerical;

        // DG Error.
        this->dg_error = std::sqrt(dot(error, dg_laplacian * error));

        // L2 Error.
        this->l2_error = std::sqrt(dot(error, mass * error));

        // L2 Errors.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {

            // Local dofs.
            std::size_t element_dofs = mesh.elements[j].dofs();

            // Global matrix indices.
            std::vector<std::size_t> indices;

            for(std::size_t k = 0; k < element_dofs; ++k)
                indices.emplace_back(j * element_dofs + k);

            this->l2_errors[j] = std::sqrt(dot(error(indices), mass(indices, indices) * error(indices)));
        }

        // Data.
        this->degree = mesh.degree;
        this->size = 0.0;
        this->elements = mesh.elements.size();

        for(const auto &element: mesh.elements)
            for(const auto &p: mesh.element(element).points)
                for(const auto &q: mesh.element(element).points) {
                    Real size = distance(p, q);
                    this->size = (size > this->size) ? size : this->size;
                }
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
        ost << "Elements: " << error.elements << std::endl;
        ost << "Degree (p): " << error.degree << std::endl;
        ost << "Size (h): " << error.size << std::endl;
        ost << "L2 Error: " << error.l2_error << std::endl;
        return ost << "DG Error: " << error.dg_error;
    }

}