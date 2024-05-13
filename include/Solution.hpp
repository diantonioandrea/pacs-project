/**
 * @file Solution.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SOLUTION_PACS
#define SOLUTION_PACS

// Type.
#include <Type.hpp>

// Mesh.
#include <Mesh.hpp>

// Vectors.
#include <Vector.hpp>

// Functor.
#include <Functor.hpp>

// Output.
#include <string>

namespace pacs {

    /**
     * @brief Readable and plottable solution.
     * 
     */
    struct Solution {

        // Vectors.
        Vector<Real> x;
        Vector<Real> y;
        Vector<Real> numerical;
        Vector<Real> exact;

        // CONSTRUCTORS.

        Solution() = delete;
        Solution(const Mesh &, const Vector<Real> &, const Functor &);

        // OUTPUT.

        void write(const std::string &);

    };

    // Solution evaluated over the quadrature nodes.

    std::array<Vector<Real>, 4> solution(const Mesh &, const Vector<Real> &, const Functor &);

    // Modal coefficients of the exact solution.

    Vector<Real> modal(const Mesh &, const Functor &);

}

#endif