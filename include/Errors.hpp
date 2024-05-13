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

// Output.
#include <iostream>

// Type.
#include <Type.hpp>

// Functor.
#include <Functor.hpp>

// Sparse matrices.
#include <Sparse.hpp>

// Mesh.
#include <Mesh.hpp>

// Containers.
#include <array>

namespace pacs {

    /**
     * @brief Error struct.
     * 
     */
    struct Error {
        
        std::size_t elements;
        std::size_t degree; // p.
        Real size; // h.

        Real dg_error;
        Real l2_error;

        // CONSTRUCTORS.

        Error(const Mesh &, const std::array<Sparse<Real>, 2> &, const Vector<Real> &, const Functor &);

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Error &);

    };

}

#endif