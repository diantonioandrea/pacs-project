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

// Sparse matrices.
#include <Sparse.hpp>

// Mesh.
#include <Mesh.hpp>

namespace pacs {

    // Laplacian matrix.

    std::array<Sparse<double>, 2> laplacian(const Mesh &);

    // Penalty coefficients.
    
    Vector<double> penalty(const Mesh &, const std::size_t &);

}

#endif