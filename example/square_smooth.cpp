/**
 * @file square_smooth.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Uniform refinement.
 * @date 2024-07-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>
#include <Laplacian.hpp>

#include "smooth.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

int main(int argc, char **argv) {

    // Degree.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " DEGREE." << std::endl;
        std::exit(-1);
    }

    std::size_t degree = static_cast<std::size_t>(std::stoi(argv[1]));

    std::ofstream output{"output/square_s_" + std::to_string(degree) + ".error"};

    output << "Square domain - uniform refinement." << "\n";

    std::cout << "Square domain - uniform refinement." << std::endl;
    std::cout << "Output under output/square_" + std::to_string(degree) + ".error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Diagrams.
    std::vector<std::vector<pacs::Polygon>> diagrams;

    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_125.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_250.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_500.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_1000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_2000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/square/square_4000.poly"));

    // Test.
    for(std::size_t j = 0; j < diagrams.size(); ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagrams[j], degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB, 1E-10);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

        // Solution structure (output).
        #ifndef NSOLUTIONS
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/square_s_" + std::to_string(degree) + "_" + std::to_string(j) + ".sol";
        solution.write(solfile);
        #endif

        // Output.
        output << "\n" << error << "\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << std::endl;
    }
}