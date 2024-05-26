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

#include <Type.hpp>
#include <Geometry.hpp>
#include <Algebra.hpp>
#include <Functor.hpp>

// Output.
#include <string>

namespace pacs {

    // Modal coefficients of the exact solution.

    Vector<Real> modal(const Mesh &, const Functor &);

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

}

#endif