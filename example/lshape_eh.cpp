/**
 * @file lshape_eh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a L-shaped domain. Element size adaptive refinement with estimator.
 * @date 2024-06-24
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

    std::ofstream output{"output/lshape_eh.error"};
    std::ofstream estimates_output{"output/lshape_eh.estimator"};

    output << "L-shaped domain - element size adaptive refinement with estimator." << "\n";
    estimates_output << "L-shaped domain - element size adaptive refinement with estimator." << "\n";

    std::cout << "L-shaped domain - element size adaptive refinement." << std::endl;
    std::cout << "Output under output/lshape_eh.error." << std::endl;

    // Domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/lshape_100.poly");

    // Polynomial degree.
    std::size_t degree = 2;

    // Refinement percentage.
    pacs::Real refine = 0.5L;

    // Mesh.
    pacs::Mesh mesh{domain, diagram, degree};

    // Test.
    for(std::size_t j = 0; j < 8; ++j) {

        // Mesh output.
        std::string polyfile = "output/lshape_eh_" + std::to_string(j) + ".poly";
        mesh.write(polyfile);

        // Matrices.
        auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

        // Forcing term.
        pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
        
        // Linear system solution.
        pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

        // Solution structure (output).
        pacs::Solution solution{mesh, numerical, exact};
        std::string solfile = "output/lshape_eh_" + std::to_string(j) + ".sol";
        solution.write(solfile);

        // Errors.
        pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact};

        // Output.
        output << "\n" << error << "\n\n";

        // Estimator.
        pacs::Estimator estimator{mesh, mass, numerical, source, dirichlet, {dirichlet_x, dirichlet_y}};

        // Output.
        estimates_output << "\n" << estimator << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";

        // Refinement.
        pacs::mesh_refine_size(mesh, estimator.estimates > refine * pacs::max(estimator.estimates));
    }
}