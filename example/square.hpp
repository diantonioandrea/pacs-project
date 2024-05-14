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
#include <Type.hpp>

// Math.
#include <cmath>
#include <cmath>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

/**
 * @brief Exact solution.
 * 
 * @return pacs::Real 
 */
pacs::Real exact(const pacs::Real &x, const pacs::Real &y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}