/**
 * @file Legendre.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>

#include <cassert>

namespace pacs {

    /**
     * @brief Binomial coefficient.
     * 
     * @param n 
     * @param m 
     * @return std::size_t 
     */
    std::size_t binomial(const std::size_t &n, const std::size_t &k) {
        #ifndef NDEBUG // Integrity check.
        assert(n >= k);
        #endif

        return ((k != 0) && (n != k)) ? binomial(n - 1, k - 1) + binomial(n - 1, k) : 1;
    }

    /**
     * @brief Returns the evaluation of a given order Legendre polynomial over x.
     * 
     * @param x 
     * @param order 
     * @return Vector<Real> 
     */
    Vector<Real> legendre(const Vector<Real> &x, const std::size_t &order) {

        // Evaluation.
        Vector<Real> evaluation{x.length, 0.0L};

        // Tabled low-order Legendre polynomials.
        switch(order) {
            case 0:
                evaluation = 1.0L;
                break;

            case 1:
                evaluation = x;
                break;

            case 2:
                evaluation = 0.5L * (3.0L * x * x - 1.0L);
                break;

            case 3:
                evaluation = 0.5L * (5.0L * x * x * x - 3.0L * x);
                break;

            case 4:
                evaluation = 0.125L * (35.0L * x * x * x * x - 30.0L * x * x + 3.0L);
                break;

            case 5:
                evaluation = 0.125L * (63.0L * x * x * x * x * x - 70.0L * x * x * x + 15.0L * x);
                break;

            case 6:
                evaluation = 0.0625L * (231.0L * x * x * x * x * x * x - 315.0L * x * x * x * x + 105.0L * x * x - 5.0L);
                break;

            default:
                break;
        }

        // Recursive formula for higher orders.
        if(order > 6) {
            for(std::size_t k = 0; k <= order; ++k) {
                Vector<Real> base = Vector<Real>(x.length, 1.0L);

                for(std::size_t h = 0; h < k; ++h)
                    base *= (x - 1.0L) / 2.0L;

                evaluation += binomial(order, k) * binomial(order + k, order) * base;
            }
        }

        return evaluation;
    }

    /**
     * @brief Returns the evaluation of the derivative of a given order Legendre polynomial over x.
     * 
     * @param x 
     * @param order 
     * @return Vector<Real> 
     */
    Vector<Real> grad_legendre(const Vector<Real> &x, const std::size_t &order) {

        // Evaluation.
        Vector<Real> evaluation{x.length, 0.0L};

        // Tabled low-order Legendre polynomials.
        switch(order) {
            case 1:
                evaluation = 1.0L;
                break;

            case 2:
                evaluation = 3.0L * x;
                break;

            case 3:
                evaluation = 0.5L * (15.0L * x * x - 3.0L);
                break;

            case 4:
                evaluation = 0.5L * (35.0L  * x * x * x - 15.0L * x);
                break;

            case 5:
                evaluation = 0.125L * (315.0L * x * x * x * x - 210.0L * x * x + 15.0L);
                break;

            case 6:
                evaluation = 0.125L * (693.0L * x * x * x * x * x - 630.0L * x * x * x + 105.0L * x);
                break;

            default:
                break;
        }

        // Recursive formula for higher orders.
        if(order > 6) {
            for(std::size_t k = 1; k <= order; ++k) {
                Vector<Real> base = Vector<Real>(x.length, 1.0L);

                for(std::size_t h = 0; h < k - 1; ++h)
                    base *= (x - 1.0L) / 2.0L;

                evaluation += k * binomial(order, k) * binomial(order + k, order) * base;
            }
        }

        return evaluation;
    }

    /**
     * @brief Returns the evaluation of the second derivative of a given order Legendre polynomial over x.
     * 
     * @param x 
     * @param order 
     * @return Vector<Real> 
     */
    Vector<Real> lap_legendre(const Vector<Real> &x, const std::size_t &order) {

        // Evaluation.
        Vector<Real> evaluation{x.length, 0.0L};

        // Tabled low-order Legendre polynomials.
        switch(order) {
            case 2:
                evaluation = 3.0L;
                break;

            case 3:
                evaluation = 15.0L * x;
                break;

            case 4:
                evaluation = 0.5L * (105.0L  * x * x - 15.0L);
                break;

            case 5:
                evaluation = 0.25L * (630.0L * x * x * x - 210.0L * x);
                break;

            case 6:
                evaluation = 0.125L * (3465.0L * x * x * x * x - 1890.0L * x * x + 105.0L);
                break;

            default:
                break;
        }

        // Recursive formula for higher orders.
        if(order > 6) {
            for(std::size_t k = 2; k <= order; ++k) {
                Vector<Real> base = Vector<Real>(x.length, 1.0L);

                for(std::size_t h = 0; h < k - 2; ++h)
                    base *= (x - 1.0L) / 2.0L;

                evaluation += k * (k - 1) * binomial(order, k) * binomial(order + k, order) * base;
            }
        }

        return evaluation;
    }

}