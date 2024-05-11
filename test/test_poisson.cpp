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

// Algebra.
#include <Algebra.hpp>

// Testing Poisson problem.
#include <Laplacian.hpp>
#include <Forcing.hpp>
#include <Solution.hpp>
using pacs::Real;

Real source(const Real &, const Real &);
Real exact(const Real &, const Real &);

int main() {

    // Constructs a mesh.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, 30)};

    // Source and exact solution.
    pacs::Functor test_source{source};
    pacs::Function test_exact{exact};

    // Writes mesh informations to a file.
    mesh.write("poisson.poly");

    // Builds the laplacian matrix.
    pacs::Sparse<Real> laplacian = pacs::laplacian(mesh)[1];

    // Builds the forcing term.
    pacs::Vector<Real> forcing = pacs::forcing(mesh, test_source);
    
    // Solution by CGM.
    pacs::Vector<Real> solution = pacs::solve(laplacian, forcing);

    // Output.
    auto [x, y, n_solution, e_solution] = pacs::solution(mesh, solution, test_exact);

}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @return Real 
 */
Real source(const Real &x, const Real &y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

/**
 * @brief Test exact solution.
 * 
 * @return Real 
 */
Real exact(const Real &x, const Real &y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}