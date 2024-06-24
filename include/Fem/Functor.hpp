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

#include <Type.hpp>
#include <Algebra.hpp>

namespace pacs {

    /**
     * @brief Zero function.
     * 
     * @return Real 
     */
    inline Real zero(const Real &, const Real &) { return 0.0; }

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
            
            Functor();
            Functor(const Function &);

            // EVALUATION.

            Real operator ()(const Real &, const Real &) const;
            Vector<Real> operator ()(const Vector<Real> &, const Vector<Real> &) const;
    };

}

#endif