/**
 * @file Functor.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FUNCTOR_PACS
#define FUNCTOR_PACS

// Type.
#include <Type.hpp>

// Vectors.
#include <Vector.hpp>

namespace pacs {

    /**
     * @brief Functor class.
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

}

#endif