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
    Error::Error(const Mesh &mesh, const std::array<Sparse<Real>, 2> &matrices, const Vector<Real> &numerical, const Functor &exact) {

        // Matrices.
        auto [mass, dg_laplacian] = matrices;

        // Error vector.
        Vector<Real> error = modal(mesh, exact) - numerical;

        // DG Error.
        this->dg_error = std::sqrt(dot(error, dg_laplacian * error));

        // L2 Error.
        this->l2_error = std::sqrt(dot(error, mass * error));

        // Data.
        this->degree = 0;
        this->size = 0.0;
        this->elements = mesh.elements.size();

        for(const auto &element: mesh.elements) {
            for(const auto &edge: element.edges)
                this->size = (std::abs(mesh.edge(edge)) > this->size) ? std::abs(mesh.edge(edge)) : this->size;

            this->degree = (element.degree > this->degree) ? element.degree : this->degree;
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