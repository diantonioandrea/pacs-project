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

    diagrams.emplace_back(pacs::mesh_diagram("data/square_1000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_2000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_3000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_4000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_5000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square_6000.poly"));

    // Polynomial degree.
    std::size_t degree = 1;

    // Test.
    for(std::size_t j = 0; j < diagrams.size(); ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagrams[j], degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB, 1E-10);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/square_" + std::to_string(j) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";
    }
}