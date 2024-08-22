/**
 * @file square.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Source and Exact solution.
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Type.
#include <Base.hpp>

// Math.
#include <cmath>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

/**
 * @brief Exact solution.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact(const pacs::Real &x, const pacs::Real &y) {
    return (1.0L - std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y) * (1.0L - x);
}

/**
 * @brief Exact solution, x derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_x(const pacs::Real &x, const pacs::Real &y) {
    return (100.0L * std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y) * (1.0L - x) - (1.0L - std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y);
}

/**
 * @brief Exact solution, y derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_y(const pacs::Real &x, const pacs::Real &y) {
    return M_PI * (1.0L - std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::cos(M_PI * y) * (1.0L - x);
}

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return (10000.0L * std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y) * (1.0L - x) + (200.0L * std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y) + M_PI * M_PI * (1.0L - std::exp(-100.0L * x)) / (1.0L - std::exp(-100.0)) * std::sin(M_PI * y) * (1.0L - x);
}