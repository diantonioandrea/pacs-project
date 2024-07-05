/**
 * @file poisson.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>
#include <Laplacian.hpp>
#include <Fem.hpp>
using namespace pacs;

// Exact solution.
inline Real exact(const Real &x, const Real &y) {
    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
}

// Source.
Real source(const Real &x, const Real &y) {
    return 8.0 * M_PI * M_PI * exact(x, y);
}

// Dirichlet BC.
inline Real dirichlet(const Real &x, const Real &y) {
    return exact(x, y);
}

int main() {

    // Domain and mesh.
    Point a{0.0, 0.0};
    Point b{1.0, 0.0};
    Point c{1.0, 1.0};
    Point d{0.0, 1.0};

    Mesh mesh{{{a, b, c, d}}, mesh_diagram("data/square/square_30.poly")};

    // Matrices.
    auto [M, A, DGA] = laplacian(mesh);

    // Forcing term.
    Vector<Real> B = forcing(mesh, source, dirichlet);

    // Numerical solution.
    Vector<Real> numerical = solve(A, B);

    // Errors.
    Error error{mesh, {M, DGA}, numerical, exact, {exact_x, exact_y}};

}