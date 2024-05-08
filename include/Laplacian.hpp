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

    Sparse<double> laplacian(const Mesh &);

}

#endif