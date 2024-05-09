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

// Testing Laplacian and Forcing.
#include <Laplacian.hpp>
#include <Forcing.hpp>

// Algebra.
#include <Algebra.hpp>

double test(const double &, const double &);

int main() {

    // Constructs a mesh.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, 32)};

    // Source.
    pacs::Source source{test};

    // Writes mesh informations to a file.
    mesh.write("poisson.poly");

    // Builds the laplacian matrix.
    pacs::Sparse<double> laplacian = pacs::laplacian(mesh)[0];

    // Builds the forcing term.
    pacs::Vector<double> forcing = pacs::forcing(mesh, source);
    
    // Solution by CGM.
    pacs::Vector<double> solution = pacs::solve(laplacian, forcing);

    // Output.
    std::cout << solution << std::endl;

}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @return double 
 */
double test(const double &x, const double &y) {
    return M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}