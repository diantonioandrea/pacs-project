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

#include <Type.hpp>
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

        // h-Adaptivity.
        Real h_estimate = 0.0L;
        Vector<Real> h_estimates;

        // p-Adaptivity.
        Real p_estimate = 0.0L;
        Vector<Real> p_estimates;

        // CONSTRUCTORS.

        Estimator(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Functor &, const Functor &dirichlet = Functor{}, const TwoFunctor &dirichlet_gradient = TwoFunctor{});

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Estimator &);

    };

}

#endif