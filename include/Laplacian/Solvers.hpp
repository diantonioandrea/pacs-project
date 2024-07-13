/**
 * @file Solvers.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LAP_SOLVERS_PACS
#define LAP_SOLVERS_PACS

#include <Type.hpp>
#include <Algebra.hpp>

#include "Laplacian.hpp"

namespace pacs {

    // Custom Laplacian solver.

    Vector<Real> lapsolver(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Real &TOL = 1E-15);

}

#endif