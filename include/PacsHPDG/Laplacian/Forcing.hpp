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

#include "../Base.hpp"
#include "../Geometry.hpp"
#include "../Fem.hpp"

namespace pacs {

    // RHS.

    Vector<Real> forcing(const Mesh &, const Functor &, const Functor &dirichlet = Functor{}, const Real &penalty_coefficient = 10.0);
}

#endif