/**
 * @file square.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Uniform refinement.
 * @date 2024-05-13
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

    std::ofstream output{"output/square.error"};

    output << "Square domain - uniform refinement." << "\n";

    std::cout << "Square domain - uniform refinement." << std::endl;
    std::cout << "Output under output/square.error." << std::endl;

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
        std::string polyfile = "output/square_" + std::to_string(j) + ".poly";
        std::vector<pacs::Polygon> diagram = std::filesystem::exists(polyfile) ? pacs::mesh_diagram(polyfile) : pacs::mesh_diagram(domain, elements);

        pacs::Mesh mesh{domain, diagram, 3};
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
        std::string surffile = "output/square_" + std::to_string(j) + ".surf";
        solution.write(surffile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << (laplacian * numerical - forcing).norm() << "\n";
    }
}