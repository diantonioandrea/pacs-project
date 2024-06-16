/**
 * @file square.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Uniform refinement.
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>
#include <Laplacian.hpp>

#include "square.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

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

    // Diagrams.
    std::vector<std::vector<pacs::Polygon>> diagrams;

    diagrams.emplace_back(pacs::mesh_diagram("data/square_100.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_200.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_400.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_800.poly"));

    // Polynomial degree.
    std::size_t degree = 2;

    // Test.
    for(std::size_t j = 0; j < diagrams.size(); ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagrams[j], degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string contfile = "output/square_" + std::to_string(j) + ".cont";
        solution.write(contfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";
    }
}