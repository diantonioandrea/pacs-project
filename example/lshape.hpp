/**
 * @file lshape.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Type.
#include <Type.hpp>

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
    pacs::Real rho2 = x * x + y * y;
    pacs::Real theta = (std::atan2(y, x) >= 0) ? std::atan2(y, x) : std::atan2(y, x) + 2.0L * M_PI;

    return std::cbrt(rho2) * std::sin(2.0 * theta / 3.0);
}

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 0.0;
}

/**
 * @brief Dirichlet boundary conditions.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real dirichlet(const pacs::Real &x, const pacs::Real &y) {
    return exact(x, y);
}