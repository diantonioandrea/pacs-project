/**
 * @file square_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement.
 * @date 2024-05-14
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

// Test functions.
#include "square.hpp"

int main() {

    std::cout << "Square domain - element size adaptive refinement." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Diagram.
    std::vector<pacs::Polygon> diagram  = pacs::mesh_diagram(domain, 16);

    // Sequence of meshes.
    for(std::size_t j = 0; j < 6; ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagram};

        // Mesh output.
        std::string polyfile = "square_h_" + std::to_string(j) + ".poly";
        mesh.write(polyfile);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = laplacian.solve<pacs::RFOM>(forcing);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string surffile = "square_h_" + std::to_string(j) + ".surf";
        solution.write(surffile);

        // Output.
        std::cout << "\n" << error << "\n" << std::endl;
        
        std::cout << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << std::endl;
        std::cout << "Residual: " << (laplacian * numerical - forcing).norm() << "\n" << std::endl;

        std::cout << "Mesh saved to " << polyfile << std::endl;
        std::cout << "Surface saved to " << surffile << std::endl;

        // Refinement.
        diagram = pacs::mesh_refine(mesh, pacs::highest(error.l2_errors, 3));
    }
}