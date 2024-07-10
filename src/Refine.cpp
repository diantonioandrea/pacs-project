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
        Real refine = 0.75L;
        Real speed = 1.5;

        // Masks.
        Mask p_mask = estimator.fits > speed;
        Mask h_mask = (estimator.estimates * estimator.estimates) > refine * sum(estimator.estimates * estimator.estimates) / mesh.elements.size();

        // Strategy.
        for(std::size_t j = 0; j < mesh.elements.size(); ++j) {
            if(!h_mask[j]) // p-Refine only error-marked elements.
                p_mask[j] = false;
                
            if(p_mask[j] && h_mask[j])
                h_mask[j] = false;
        }

        // Refinements.
        mesh_refine_degree(mesh, p_mask);
        mesh_refine_size(mesh, h_mask);
    }

}