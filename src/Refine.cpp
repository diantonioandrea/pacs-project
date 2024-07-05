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
        Real refine = 0.75;
        Real speed = 1.0;

        // Estimates, maximum estimate and fits.
        Vector<Real> estimates = estimator.estimates;
        Real maximum = max(estimates);
        Vector<Real> fits = estimator.fits;

        // Mask.
        Mask mask = estimates > refine * maximum;

        // Masks.
        Mask p_mask(mask.size(), false);
        Mask h_mask(mask.size(), false);

        for(std::size_t j = 0; j < mask.size(); ++j)
            if(mask[j]) {
                if(fits[j] > speed)
                    p_mask[j] = true;
                else
                    h_mask[j] = true;
            }

        // Refinements.
        mesh_refine_degree(mesh, p_mask);
        mesh_refine_size(mesh, h_mask);

    }

}