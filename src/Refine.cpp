/**
 * @file Refine.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-27
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Laplacian.hpp>

namespace pacs {

    /**
     * @brief hp-Adaptively refine a mesh.
     * 
     * @param mesh 
     * @param estimator 
     */
    void mesh_refine(Mesh &mesh, const Estimator &estimator) {

        // Refinement steps.
        std::array<Real, 2> steps{0.5, 0.75};

        // Estimates and maximum estimate.
        Vector<Real> estimates = estimator.estimates;
        Real maximum = max(estimates);

        // p-Mask.
        Mask p_mask = Mask(estimates.length);
        Mask lower_p_mask = estimates > steps[0] * maximum;
        Mask upper_p_mask = estimates < steps[1] * maximum;

        for(std::size_t j = 0; j < estimates.length; ++j)
            p_mask[j] = lower_p_mask[j] && upper_p_mask[j];

        // h-Mask.
        Mask h_mask = estimates > steps[1] * maximum;

        // Refinements.
        mesh_refine_degree(mesh, p_mask);
        mesh_refine_size(mesh, h_mask);

    }

}