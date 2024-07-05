/**
 * @file test_poisson.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Laplacian.hpp>
#include <Fem.hpp>

#include <iostream>

pacs::Real exact(const pacs::Real &, const pacs::Real &);
pacs::Real exact_x(const pacs::Real &, const pacs::Real &);
pacs::Real exact_y(const pacs::Real &, const pacs::Real &);
pacs::Real source(const pacs::Real &, const pacs::Real &);
pacs::Real dirichlet(const pacs::Real &, const pacs::Real &);

int main() {

    // Loads a mesh.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram("data/square/square_30.poly"), 3};

    // Writes mesh informations to a file.
    mesh.write("output/poisson.poly");

    // Builds the laplacian matrix.
    auto [mass, laplacian, dg_laplacian] = pacs::laplacian(mesh);

    // Builds the forcing term.
    pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);

    // Linear system solution.
    pacs::Vector<pacs::Real> numerical = pacs::solve(laplacian, forcing, pacs::BICGSTAB);

    // Errors.
    pacs::Error error{mesh, {mass, dg_laplacian}, numerical, exact, {exact_x, exact_y}};

    // Solution structure (output).
    pacs::Solution solution{mesh, numerical, exact};
    solution.write("output/poisson.sol");

    // Output.
    std::cout << "\n" << error << "\n" << std::endl;

}

/**
 * @brief Exact solution.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact(const pacs::Real &x, const pacs::Real &y) {
    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
}

/**
 * @brief Exact solution, x derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_x(const pacs::Real &x, const pacs::Real &y) {
    return 2.0L * M_PI * std::cos(2.0L * M_PI * x) * std::cos(2.0L * M_PI * y);
}

/**
 * @brief Exact solution, y derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_y(const pacs::Real &x, const pacs::Real &y) {
    return -2.0L * M_PI * std::sin(2.0L * M_PI * x) * std::sin(2.0L * M_PI * y);
}

/**
 * @brief Test source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 8.0 * M_PI * M_PI * exact(x, y);
}

/**
 * @brief Dirichlet BC.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real dirichlet(const pacs::Real &x, const pacs::Real &y) {
    return exact(x, y);
}