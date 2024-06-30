/**
 * @file square_eh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement with estimator.
 * @date 2024-06-24
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

int main(int argc, char **argv) {

    // Degree.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " DEGREE." << std::endl;
        std::exit(-1);
    }

    std::size_t degree = static_cast<std::size_t>(std::stoi(argv[1]));

    std::ofstream output{"output/square_eh_" + std::to_string(degree) + ".error"};
    std::ofstream estimates_output{"output/square_eh_" + std::to_string(degree) + ".estimator"};

    output << "Square domain - element size adaptive refinement with estimator." << "\n";
    estimates_output << "Square domain - element size adaptive refinement with estimator." << "\n";

    std::cout << "Square domain - element size adaptive refinement with estimator." << std::endl;
    std::cout << "Output under output/square_eh_DEGREE.error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/square_100.poly");

    // Refinement percentage.
    pacs::Real refine = 0.25L;

    // Mesh.
    pacs::Mesh mesh{domain, diagram, degree};

    // Sequence of meshes.
    for(std::size_t index = 0; index < 10; ++index) {

        // Verbosity.
        std::cout << "\nDEGREE: " << degree << "\nINDEX: " << index << "\n" << std::endl;

        // Mesh output.
        std::string polyfile = "output/square_eh_" + std::to_string(degree) + "_" + std::to_string(index) + ".poly";
        mesh.write(polyfile);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // // Solution structure (output).
        // pacs::Solution solution{mesh, numerical, exact};
        // std::string solfile = "output/square_eh_" + std::to_string(degree) + "_" + std::to_string(index) + ".sol";
        // solution.write(solfile);

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

        // Exit.
        if(index == 9)
            break;

        // Refinement.
        pacs::mesh_refine_size(mesh, estimator.estimates > refine * pacs::max(estimator.estimates));
    }
}
