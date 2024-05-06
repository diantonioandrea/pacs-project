/**
 * @file Quadrature.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef QUADRATURE_PACS
#define QUADRATURE_PACS

// Vectors.
#include <Vector.hpp>

// Geometry.
#include <Geometry.hpp>

// Containers.
#include <vector>

// Quadrature tolerance.
#ifndef QUADRATURE_TOLERANCE
#define QUADRATURE_TOLERANCE 1E-12
#endif

namespace pacs {

    // Gauss-Legendre quadrature nodes.

    std::vector<Vector<double>> gauss_legendre(const double &, const double &, const std::size_t &);

}

#endif