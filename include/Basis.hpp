/**
 * @file Basis.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASIS_PACS
#define BASIS_PACS

// Type.
#include <Type.hpp>

// Vectors.
#include <Vector.hpp>

// Matrices.
#include <Matrix.hpp>

// Mesh.
#include <Mesh.hpp>

// Containers.
#include <array>

namespace pacs {

    // Basis functions.

    std::array<Matrix<Real>, 3> basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);

}

#endif