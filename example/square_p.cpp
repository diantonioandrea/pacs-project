/**
 * @file square_p.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Uniform degree refinement.
 * @date 2024-06-03
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

    std::ofstream output{"output/square_p.error"};

    output << "Square domain - uniform degree refinement." << "\n";

    std::cout << "Square domain - uniform degree refinement." << std::endl;
    std::cout << "Output under output/square_p.error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/square_1000.poly");

    // Test.
    for(std::size_t degree = 1; degree <= 6; ++degree) {

        // Mesh.
        pacs::Mesh mesh{domain, diagram, degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/square_p_" + std::to_string(degree) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";
    }
}