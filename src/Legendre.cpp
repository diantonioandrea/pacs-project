/**
 * @file Legendre.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Legendre.hpp>

// Assertions.
#include <cassert>

namespace pacs {

    /**
     * @brief Returns the evaluation of a given order Legendre polynomials over x.
     * 
     * @param x 
     * @param order 
     * @return Vector<Real> 
     */
    Vector<Real> legendre(const Vector<Real> &x, const std::size_t &order) {
        #ifndef NDEBUG // Testing.
        assert(order < 7);
        #endif

        // Evaluation.
        Vector<Real> evaluation{x.length, 0.0};

        // Tabled low-order Legendre polynomials.
        switch(order) {
            case 0:
                evaluation = 1.0;
                break;

            case 1:
                evaluation = x;
                break;

            case 2:
                evaluation = 0.5 * (3.0 * x * x - 1.0);
                break;

            case 3:
                evaluation = 0.5 * (5.0 * x * x * x - 3.0 * x);
                break;

            case 4:
                evaluation = 0.125 * (35.0 * x * x * x * x - 30.0 * x * x + 3.0);
                break;

            case 5:
                evaluation = 0.125 * (63.0 * x * x * x * x * x - 70.0 * x * x * x + 15.0 * x);
                break;

            case 6:
                evaluation = 0.125 * (231.0 * x * x * x * x * x * x - 315.0 * x * x * x * x + 105.0 * x * x - 5.0);
                break;

            default:
                break;
        }

        return evaluation;
    }

    /**
     * @brief Returns the evaluation of the derivative of a given order Legendre polynomials over x.
     * 
     * @param x 
     * @param order 
     * @return Vector<Real> 
     */
    Vector<Real> grad_legendre(const Vector<Real> &x, const std::size_t &order) {
        #ifndef NDEBUG // Testing.
        assert(order < 7);
        #endif

        // Evaluation.
        Vector<Real> evaluation{x.length, 0.0};

        // Tabled low-order Legendre polynomials.
        switch(order) {
            case 1:
                evaluation = 1.0;
                break;

            case 2:
                evaluation = 3.0 * x;
                break;

            case 3:
                evaluation = 0.5 * (15.0 * x * x - 3.0);
                break;

            case 4:
                evaluation = 0.5 * (35.0  * x * x * x - 15.0 * x);
                break;

            case 5:
                evaluation = 0.125 * (315.0 * x * x * x * x - 210.0 * x * x + 15.0);
                break;

            case 6:
                evaluation = 0.250 * (693.0 * x * x * x * x * x - 630.0 * x * x * x + 105.0 * x);
                break;

            default:
                break;
        }

        return evaluation;
    }

}