/**
 * @file Solution.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SOLUTION_PACS
#define SOLUTION_PACS

// Type.
#include <Type.hpp>

// Mesh.
#include <Mesh.hpp>

// Vectors.
#include <Vector.hpp>

// Functor.
#include <Functor.hpp>

namespace pacs {

    // Solution evaluated over the quadrature nodes.

    std::array<Vector<Real>, 4> solution(const Mesh &, const Vector<Real> &, const Functor &);

}

#endif