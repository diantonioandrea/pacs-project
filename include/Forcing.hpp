/**
 * @file Forcing.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FORCING_PACS
#define FORCING_PACS

// Type.
#include <Type.hpp>

// Mesh.
#include <Mesh.hpp>

// Functor.
#include <Functor.hpp>

namespace pacs {

    // RHS.

    Vector<Real> forcing(const Mesh &, const Functor &);
}

#endif