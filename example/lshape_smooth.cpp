/**
 * @file lshape_smooth.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a L-shaped domain. Uniform refinement.
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

    std::ofstream output{"output/lshape_s_" + std::to_string(degree) + ".error"};

    output << "L-shaped domain - uniform refinement." << "\n";

    std::cout << "L-shaped domain - uniform refinement." << std::endl;
    std::cout << "Output under output/lshape_" + std::to_string(degree) + ".error." << std::endl;

    // Domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};

    // Diagrams.
    std::vector<std::vector<pacs::Polygon>> diagrams;

    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_125.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_250.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_500.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_1000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_2000.poly"));
    diagrams.emplace_back(pacs::mesh_diagram("data/lshape/lshape_4000.poly"));

    // Test.
    for(std::size_t j = 0; j < diagrams.size(); ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagrams[j], degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

        // Solution structure (output).
        #ifndef NSOLUTIONS
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/lshape_s_" + std::to_string(degree) + "_" + std::to_string(j) + ".sol";
        solution.write(solfile);
        #endif

        // Output.
        output << "\n" << error << "\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << std::endl;
    }
}