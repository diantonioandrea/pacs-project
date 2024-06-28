/**
 * @file square_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement.
 * @date 2024-05-31
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

    std::ofstream output{"output/square_h.error"};

    output << "Square domain - element size adaptive refinement." << "\n";

    std::cout << "Square domain - element size adaptive refinement." << std::endl;
    std::cout << "Output under output/square_h.error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/square_1000.poly");

    // Polynomial degree.
    std::size_t degree = 2;

    // Refinement percentage.
    pacs::Real refine = 0.25L;

    // Mesh.
    pacs::Mesh mesh{domain, diagram, degree};

    // Sequence of meshes.
    while(mesh.elements.size() < 20000) {

        // Mesh output.
        std::string polyfile = "output/square_h_" + std::to_string(mesh.elements.size()) + ".poly";
        mesh.write(polyfile);

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
        std::string solfile = "output/square_h_" + std::to_string(j) + ".sol";
        solution.write(solfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << pacs::norm(laplacian * numerical - forcing) << "\n";

        // Refinement.
        pacs::mesh_refine_size(mesh, error.l2_errors > refine * pacs::max(error.l2_errors));
    }
}
