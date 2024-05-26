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

#include <Type.hpp>

#include <Algebra.hpp>
#include <Geometry.hpp>

#include <vector>

// Quadrature tolerance.
#ifndef QUADRATURE_TOLERANCE
#define QUADRATURE_TOLERANCE 1E-14
#endif

namespace pacs {

    // Gauss-Legendre quadrature nodes.

    std::array<Vector<Real>, 2> gauss_legendre(const Real &, const Real &, const std::size_t &);

    // Reference interval and square quadrature nodes.

    std::array<Vector<Real>, 2> quadrature_1d(const std::size_t &);
    std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &);

}

#endif