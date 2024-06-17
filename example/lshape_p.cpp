/**
 * @file lshape_p.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a L-shaped domain. Uniform degree refinement.
 * @date 2024-06-17
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>
#include <Laplacian.hpp>

#include "lshape.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

int main() {

    std::ofstream output{"output/lshape_p.error"};

    output << "L-shaped domain - uniform degree refinement." << "\n";

    std::cout << "L-shaped domain - uniform degree refinement." << std::endl;
    std::cout << "Output under output/lshape_p.error." << std::endl;

    // Domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};

    // Diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/lshape_100.poly");

    // Test.
    for(std::size_t degree = 1; degree <= 6; ++degree) {

        // Mesh.
        pacs::Mesh mesh{domain, diagram, degree};

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/lshape_p_" + std::to_string(degree) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";
    }
}