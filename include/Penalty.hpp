/**
 * @file Penalty.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PENALTY_PACS
#define PENALTY_PACS

// Type.
#include <Type.hpp>

// Vectors.
#include <Algebra.hpp>

// Mesh.
#include <Mesh.hpp>

namespace pacs {

    // Penalty coefficients.
    
    Vector<Real> penalty(const Mesh &, const std::size_t &);

}

#endif