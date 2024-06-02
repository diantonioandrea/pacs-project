/**
 * @file square_comparison_h.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Element size adaptive refinement and comparison.
 * @date 2024-05-31
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#define ALGEBRA_TOLERANCE 1E-8 // Faster.

#include <Fem.hpp>
#include <Laplacian.hpp>

#include "square.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

int main() {

    std::ofstream output{"output/square.error"}, output_h{"output/square_h.error"};

    output << "Square domain - element size adaptive refinement and comparison." << "\n";
    output_h << "Square domain - element size adaptive refinement and comparison." << "\n";

    std::cout << "Square domain - element size adaptive refinement and comparison." << std::endl;
    std::cout << "Output under output/square.error and output/square_h.error." << std::endl;

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Initial diagram.
    std::vector<pacs::Polygon> diagram = pacs::mesh_diagram("data/square_200.poly");

    // Elements.
    std::vector<std::size_t> elements;

    // Polynomial degree.
    std::size_t degree = 3;

    // Refinement percentage.
    pacs::Real refine = 0.5L;

    // Sequence of meshes.
    for(std::size_t j = 0; j < 4; ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, diagram, degree};
        elements.emplace_back(mesh.elements.size());

        // Mesh output.
        std::string polyfile = "output/square_h_" + std::to_string(j) + ".poly";
        mesh.write(polyfile);

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
        std::string contfile = "output/square_h_" + std::to_string(j) + ".cont";
        solution.write(contfile);

        // Output.
        output_h << "\n" << error << "\n\n";
        
        output_h << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output_h << "Residual: " << norm(laplacian * numerical - forcing) << "\n";

        // Refinement.
        diagram = pacs::mesh_refine(mesh, error.l2_errors > refine * pacs::max(error.l2_errors));
    }

    // Closes output_h.
    output_h.close();

    // Comparison.

    // Sequence of meshes.
    for(std::size_t j = 0; j < elements.size(); ++j) {

        // Mesh.
        pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, elements[j]), degree};

        // Mesh output.
        std::string polyfile = "output/square_" + std::to_string(j) + ".poly";
        mesh.write(polyfile);

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
        std::string contfile = "output/square_" + std::to_string(j) + ".cont";
        solution.write(contfile);

        // Output.
        output << "\n" << error << "\n\n";
        
        output << "Laplacian: " << laplacian.rows << " x " << laplacian.columns << "\n";
        output << "Residual: " << norm(laplacian * numerical - forcing) << "\n";
    }
}
