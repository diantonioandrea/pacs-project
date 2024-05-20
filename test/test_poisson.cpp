/**
 * @file test_poisson.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Testing Poisson problem.
#include <Laplacian.hpp>
#include <Forcing.hpp>
#include <Solution.hpp>

pacs::Real source(const pacs::Real &, const pacs::Real &);

int main() {

    // Loads a mesh.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram("data/example.poly")};

    // Source and exact solution.
    pacs::Functor test_source{source};

    // Writes mesh informations to a file.
    mesh.write("poisson.poly");

    // Builds the laplacian matrix.
    auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

    // Builds the forcing term.
    pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, test_source);

    // Linear system solution.
    pacs::Vector<pacs::Real> solution = laplacian.solve<pacs::RFOM>(forcing);

}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}