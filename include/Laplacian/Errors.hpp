/**
 * @file Errors.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ERRORS_PACS
#define ERRORS_PACS

#include <Type.hpp>
#include <Fem.hpp>
#include <Algebra.hpp>
#include <Geometry.hpp>

#include <iostream>
#include <array>

namespace pacs {

    /**
     * @brief Error struct.
     * 
     */
    struct Error {
        
        std::size_t elements;
        std::size_t dofs;
        std::size_t degree; // p.
        Real size; // h.

        Real dg_error;
        Real l2_error;

        Vector<Real> l2_errors;
        Vector<Real> h1_errors;

        // CONSTRUCTORS.

        Error(const Mesh &, const std::array<Sparse<Real>, 2> &, const Vector<Real> &, const Functor &, const TwoFunctor &);

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Error &);

    };

}

#endif