/**
 * @file Functor.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Pacs.hpp>

namespace pacs {

    // FUNCTOR.

    // CONSTRUCTORS.

    /**
     * @brief Construct a new Functor with the zero function.
     * 
     */
    Functor::Functor(): function{zero} {}

    /**
     * @brief Construct a new Functor from a given Function.
     * 
     * @param function 
     */
    Functor::Functor(const Function &function): function{function} {}

    // EVALUATION.
    
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

    // TWOFUNCTOR.

    // CONSTRUCTORS.

    /**
     * @brief Construct a new TwoFunctor with two zero functions.
     * 
     */
    TwoFunctor::TwoFunctor(): first{zero}, second{zero} {}

    /**
     * @brief Construct a new TwoFunctor from two given Function.
     * 
     * @param first 
     * @param second 
     */
    TwoFunctor::TwoFunctor(const Function &first, const Function &second): first{first}, second{second} {}

    // EVALUATION.

    /**
     * @brief Evaluates the functions over a given point.
     * 
     * @param x 
     * @param y 
     * @return Vector<Real> 
     */
    std::array<Real, 2> TwoFunctor::operator ()(const Real &x, const Real &y) const {
        return {first(x, y), second(x, y)};
    }

    /**
     * @brief Evaluates the function over a given vector of points.
     * 
     * @param x 
     * @param y 
     * @return std::array<Vector<Real>, 2> 
     */
    std::array<Vector<Real>, 2> TwoFunctor::operator ()(const Vector<Real> &x, const Vector<Real> &y) const {
        #ifndef NDEBUG // Integrity check.
        assert(x.length == y.length);
        #endif

        Vector<Real> output_first{x.length};
        Vector<Real> output_second{x.length};

        for(std::size_t j = 0; j < x.length; ++j) {
            output_first[j] = this->first(x[j], y[j]);
            output_second[j] = this->second(x[j], y[j]);
        }

        return {output_first, output_second};
    }

}