/**
 * @file lshape_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a L-shaped domain. Element size adaptive refinement.
 * @date 2024-06-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include "lshape.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

int main(int argc, char **argv) {

    // Degree.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " DEGREE [ELEMENTS]." << std::endl;
        std::exit(-1);
    }

    std::size_t degree = static_cast<std::size_t>(std::stoi(argv[1]));
    
    // Initial diagram.
    std::size_t elements = 125;

    if(argc == 3)
        elements = static_cast<std::size_t>(std::stoi(argv[2]));

    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/lshape/lshape_" + std::to_string(elements) + ".poly");

    // "Splash".
    std::ofstream output{"output/lshape_h_" + std::to_string(elements) + "@" + std::to_string(degree) + ".error"};

    output << "L-shaped domain - element size adaptive refinement." << "\n";

    std::cout << "L-shaped domain - element size adaptive refinement." << std::endl;
    std::cout << "Output under output/lshape_h_" + std::to_string(degree) + ".error." << std::endl;

    // Domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};

    // Refinement percentage.
    pacs::Real refine = 0.75L;

    // Mesh.
    pacs::Mesh mesh{domain, diagram, degree};

    // Test.
    for(std::size_t index = 0; index < TESTS_MAX; ++index) {

        // Verbosity.
        std::cout << "\nDEGREE: " << degree << "\nINDEX: " << index << "\n" << std::endl;

        // Mesh output.
        std::string polyfile = "output/lshape_h_" + std::to_string(elements) + "@" + std::to_string(degree) + "_" + std::to_string(index) + ".poly";
        mesh.write(polyfile);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::lapsolver(mesh, laplacian, forcing);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

        // Solution structure (output).
        #ifndef NSOLUTIONS
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/lshape_h_" + std::to_string(elements) + "@" + std::to_string(degree) + "_" + std::to_string(index) + ".sol";
        solution.write(solfile);
        #endif

        // Output.
        output << "\n" << error << "\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << std::endl;

        // Exit.
        if(error.dofs > DOFS_MAX)
            break;

        // Refinement.
        pacs::mesh_refine_size(mesh, error.l2_errors > refine * pacs::max(error.l2_errors));
    }
}