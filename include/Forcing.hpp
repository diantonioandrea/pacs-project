/**
 * @file Forcing.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FORCING_PACS
#define FORCING_PACS

// Type.
#include <Type.hpp>

// Mesh.
#include <Mesh.hpp>

namespace pacs {

    /**
     * @brief Functor class for the Poisson problem.
     * 
     */
    class Functor {
        private:
            
            // Function.
            Function function;

        public:

            // CONSTRUCTORS.

            Functor(const Function &);

            // EVALUATION.

            Real operator ()(const Real &, const Real &) const;
            Vector<Real> operator ()(const Vector<Real> &, const Vector<Real> &) const;
    };

    // RHS.

    Vector<Real> forcing(const Mesh &, const Functor &);
}

#endif