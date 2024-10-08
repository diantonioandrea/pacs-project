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

#include "../Base.hpp"
#include "../Algebra.hpp"

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

    /**
     * @brief TwoFunctor class.
     * 
     */
    class TwoFunctor {
        private:

            // Functions.
            Function first;
            Function second;

        public:

            // CONSTRUCTORS.

            TwoFunctor();
            TwoFunctor(const Function &, const Function &);

            // EVALUATION.

            std::array<Real, 2> operator ()(const Real &, const Real &) const ;
            std::array<Vector<Real>, 2> operator ()(const Vector<Real> &, const Vector<Real> &) const;

    };

}

#endif