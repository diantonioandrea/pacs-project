/**
 * @file square_hp.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. hp-Adaptive refinement with estimator.
 * @date 2024-06-27
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

    std::ofstream output{"output/square_hp.error"};
    std::ofstream estimates_output{"output/square_hp.estimator"};

    output << "Square domain - hp-adaptive refinement with estimator." << "\n";
    estimates_output << "Square domain - hp-adaptive refinement with estimator." << "\n";

    std::cout << "Square domain - hp-adaptive refinement with estimator." << std::endl;
    std::cout << "Output under output/square_hp.estimator." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/square_100.poly");

    // Mesh.
    pacs::Mesh mesh{domain, diagram, 2};

    // Sequence of meshes.
    for(std::size_t j = 0; j < 6; ++j) {

        // Mesh output.
        std::string polyfile = "output/square_hp_" + std::to_string(j) + ".poly";
        mesh.write(polyfile, true);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/square_hp_" + std::to_string(j) + ".sol";
        solution.write(solfile);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Output.
        output << "\n" << error << "\n\n";

        // Estimator.
        pacs::Estimator estimator{mesh, mass, numerical, source};

        // Output.
        estimates_output << "\n" << estimator << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";

        // Refinement.
        pacs::mesh_refine(mesh, estimator);
    }
}
