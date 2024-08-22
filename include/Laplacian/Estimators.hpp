/**
 * @file Estimators.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ESTIMATORS_PACS
#define ESTIMATORS_PACS

#include <Base.hpp>
#include <Fem.hpp>
#include <Algebra.hpp>
#include <Geometry.hpp>

#include <iostream>
#include <array>

namespace pacs {

    /**
     * @brief Error estimator struct.
     * 
     */
    struct Estimator {

        // DOFs.
        std::size_t dofs;

        // Estimates.
        Real estimate = 0.0L;
        Vector<Real> estimates;

        // Fits.
        Vector<Real> fits;

        // CONSTRUCTORS.

        Estimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Functor &, const Functor &dirichlet = Functor{}, const TwoFunctor &dirichlet_gradient = TwoFunctor{});

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Estimator &);

    };

}

#endif