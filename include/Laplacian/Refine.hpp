/**
 * @file Refine.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-27
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef REFINE_PACS
#define REFINE_PACS

#include <Base.hpp>
#include <Geometry.hpp>

#include "Estimators.hpp"

namespace pacs {

    // Adaptive refinement.

    void mesh_refine(Mesh &, const Estimator &, const Real &refine = 0.75, const Real &speed = 1.0);

}

#endif