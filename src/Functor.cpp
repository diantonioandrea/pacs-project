/**
 * @file Functor.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Functor.hpp>

namespace pacs {

    // CONSTRUCTORS

    /**
     * @brief Construct a new Functor from a given Function.
     * 
     * @param function 
     */
    Functor::Functor(const Function &function): function{function} {}

    // EVALUATION
    
    /**
     * @brief Evaluates the function over a given point.
     * 
     * @param x 
     * @param y 
     * @return Real 
     */
    Real Functor::operator ()(const Real &x, const Real &y) const {
        return this->function(x, y);
    }

    /**
     * @brief Evaluates the function over a given vector of points.
     * 
     * @param x 
     * @param y 
     * @return Vector<Real> 
     */
    Vector<Real> Functor::operator ()(const Vector<Real> &x, const Vector<Real> &y) const {
        #ifndef NDEBUG // Integrity check.
        assert(x.length == y.length);
        #endif

        Vector<Real> output{x.length};

        for(std::size_t j = 0; j < output.length; ++j)
            output[j] = this->function(x[j], y[j]);

        return output;
    }

}