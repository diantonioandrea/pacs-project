/**
 * @file square.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain.
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Needed headers.
#include <Laplacian.hpp>
#include <Forcing.hpp>
#include <Solution.hpp>

pacs::Real source(const pacs::Real &, const pacs::Real &);
pacs::Real exact(const pacs::Real &, const pacs::Real &);

int main() {

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Sequence of meshes.
    for(std::size_t j = 0; j < 4; ++j) {
        std::size_t elements = std::pow(2, 4 + j);

        // Mesh.
        pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, elements)};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = laplacian.solve(forcing);

        // Error vector.
        pacs::Vector<pacs::Real> error = pacs::modal(mesh, exact) - numerical;

        // DG Error.
        pacs::Real dg_error = std::sqrt(pacs::dot(error, dg_laplacian * error));

        // L2 Error.
        pacs::Real l2_error = std::sqrt(pacs::dot(error, mass * error));

        // Output.
        std::cout << ((j != 0) ? "\n" : "") << "Elements: " << elements << std::endl;
        std::cout << "L2 Error: " << l2_error << std::endl;
        std::cout << "DG Error: " << dg_error << std::endl;
    }
}

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

/**
 * @brief Exact solution.
 * 
 * @return pacs::Real 
 */
pacs::Real exact(const pacs::Real &x, const pacs::Real &y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}