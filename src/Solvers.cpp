/**
 * @file Solvers.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Pacs.hpp>

namespace pacs {

    /**
     * @brief Custom Laplacian solver.
     * 
     * @param mesh 
     * @param A 
     * @param b 
     * @param TOL 
     * @return Vector<Real> 
     */
    Vector<Real> lapsolver(const Mesh &mesh, const Sparse<Real> &A, const Vector<Real> &b, const Real &TOL) {
        // Mass blocks.
        auto blocks = block_mass(mesh);
        
        // Solves using BICGSTAB and DBI preconditioner.
        return solve(A, b, blocks, GMRES, DBI, TOL);
    }

}