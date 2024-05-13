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

// Matrices and RHS.
#include <Laplacian.hpp>
#include <Forcing.hpp>

// Error evaluation.
#include <Errors.hpp>

// Solution plot.
#include <Solution.hpp>

// Filename.
#include <string>

pacs::Real source(const pacs::Real &, const pacs::Real &);
pacs::Real exact(const pacs::Real &, const pacs::Real &);

int main() {

    std::cout << "Square domain - uniform refinement." << std::endl;

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

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string filename = "square_" + std::to_string(mesh.elements.size()) + ".surf";
        solution.write(filename);

        // Output.
        std::cout << "\n" << error << std::endl;
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