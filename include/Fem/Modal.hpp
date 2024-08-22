/**
 * @file Modal.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MODAL_PACS
#define MODAL_PACS

#include <Base.hpp>
#include <Geometry.hpp>
#include <Algebra.hpp>

#include "Functor.hpp"

#include <string>

namespace pacs {

    // Modal coefficients of a function.

    Vector<Real> modal(const Mesh &, const Functor &);

}

#endif