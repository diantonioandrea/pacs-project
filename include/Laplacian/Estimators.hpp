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

        std::size_t elements;
        std::size_t dofs;

        // Estimates.
        Real estimate;
        Vector<Real> estimates;

        // CONSTRUCTORS.

        Estimator(const Mesh &, const std::array<Sparse<Real>, 2> &, const Vector<Real> &, const Functor &, const Functor &dirichlet = Functor{});

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Estimator &);

    };

}

#endif