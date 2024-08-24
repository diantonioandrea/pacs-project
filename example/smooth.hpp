/**
 * @file square_smooth.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a square domain. Smooth Source and Exact solution.
 * @date 2024-07-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Pacs.hpp>

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
    return std::sin(2.0L * M_PI * x) * std::cos(2.0L * M_PI * y);
}

/**
 * @brief Exact solution, x derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_x(const pacs::Real &x, const pacs::Real &y) {
    return 2.0L * M_PI * std::cos(2.0L * M_PI * x) * std::cos(2.0L * M_PI * y);
}

/**
 * @brief Exact solution, y derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_y(const pacs::Real &x, const pacs::Real &y) {
    return -2.0L * M_PI * std::sin(2.0L * M_PI * x) * std::sin(2.0L * M_PI * y);
}

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 8.0L * M_PI * M_PI * exact(x, y);
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