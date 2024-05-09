/**
 * @file Laplacian.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LAPLACIAN_PACS
#define LAPLACIAN_PACS

// Type.
#include <Type.hpp>

// Sparse matrices.
#include <Sparse.hpp>

// Mesh.
#include <Mesh.hpp>

namespace pacs {

    // Laplacian matrix.

    std::array<Sparse<Real>, 2> laplacian(const Mesh &);

    // Penalty coefficients.
    
    Vector<Real> penalty(const Mesh &, const std::size_t &);

}

#endif