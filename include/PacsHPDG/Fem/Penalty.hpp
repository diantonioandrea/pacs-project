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

#include "../Base.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

namespace pacs {

    // Penalty coefficients.
    
    Vector<Real> penalty(const Mesh &, const std::size_t &, const Real &);

}

#endif