/**
 * @file hp_refine.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-10
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>
#include <Laplacian.hpp>
#include <Fem.hpp>
using namespace pacs;

// Exact solution and source.
#include "../example/square.hpp"

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
    Vector<Real> B = forcing(mesh, source);

    // Numerical solution.
    Vector<Real> numerical = lapsolver(mesh, A, B);

    // Estimates.
    Estimator est{mesh, M, numerical, source};

    // Refinement.
    mesh_refine(mesh, est);
}