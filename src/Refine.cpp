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
     * @param refine Refinement percentage.
     * @param speed Solution's smoothness.
     */
    void mesh_refine(Mesh &mesh, const Estimator &estimator, const Real &refine, const Real &speed) {
        #ifndef NDEBUG // Integrity check.
        assert((refine > 0.0L) && (refine < 1.0L));
        assert(speed > 0.0L);
        #endif

        // Masks.
        Mask p_mask = estimator.fits > speed;
        Mask h_mask = (estimator.estimates * estimator.estimates) > refine * sum(estimator.estimates * estimator.estimates) / mesh.elements.size();

        // Strategy.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {
            if(!h_mask[j]) // p-Refine only error-marked elements.
                p_mask[j] = false;
                
            if(p_mask[j] && h_mask[j]) // p > h.
                h_mask[j] = false;
        }

        // Refinements.
        mesh_refine_degree(mesh, p_mask);
        mesh_refine_size(mesh, h_mask);
    }

}