/**
 * @file Laplacian.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LAPLACIAN_MATRIX__PACS
#define LAPLACIAN_MATRIX__PACS

#include <Type.hpp>
#include <Algebra.hpp>
#include <Geometry.hpp>

namespace pacs {

    // Laplacian matrix.

    std::array<Sparse<Real>, 3> laplacian(const Mesh &);

    // Blocks.

    std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &);
}

#endif