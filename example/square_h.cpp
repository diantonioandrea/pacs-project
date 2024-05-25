/**
 * @file square_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement.
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// IO handling.
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

// Matrices and RHS.
#include <Laplacian.hpp>
#include <Forcing.hpp>

// Error evaluation.
#include <Errors.hpp>

// Solution plot.
#include <Solution.hpp>

// Test functions.
#include "square.hpp"

int main() {

    std::ofstream output{"output/square_h.error"};

    output << "Square domain - element size adaptive refinement." << "\n";

    std::cout << "Square domain - element size adaptive refinement." << std::endl;
    std::cout << "Output under output/square_h.error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 16);

    // Sequence of meshes.
    for(std::size_t j = 0; j < 6; ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagram, 3};

        // Mesh output.
        std::string polyfile = "output/square_h_" + std::to_string(j) + ".poly";
        mesh.write(polyfile);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = laplacian.solve<pacs::GMRES>(forcing);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string contfile = "output/square_h_" + std::to_string(j) + ".cont";
        solution.write(contfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << (laplacian * numerical - forcing).norm() << "\n";

        // Refinement.
        diagram = pacs::mesh_refine(mesh, pacs::highest(error.l2_errors, 3));
    }
}
