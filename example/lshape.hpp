/**
 * @file lshape.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Poisson on a L-shaped domain. Source and Exact solution.
 * @date 2024-06-16
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
    pacs::Real rho = std::sqrt(x * x + y * y);
    pacs::Real theta = (std::atan2(y, x) >= 0) ? std::atan2(y, x) : std::atan2(y, x) + 2.0L * M_PI;

    return std::pow(rho, 2.0L / 3.0L) * std::sin(2.0 * theta / 3.0);
}

/**
 * @brief Exact solution, x derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_x(const pacs::Real &x, const pacs::Real &y) {
    pacs::Real rho = std::sqrt(x * x + y * y);
    pacs::Real theta = (std::atan2(y, x) >= 0) ? std::atan2(y, x) : std::atan2(y, x) + 2.0L * M_PI;

    return -(2.0L / 3.0L) * std::pow(rho, -1.0L / 3.0L) * std::sin(theta / 3.0L);
}

/**
 * @brief Exact solution, y derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real exact_y(const pacs::Real &x, const pacs::Real &y) {
    pacs::Real rho = std::sqrt(x * x + y * y);
    pacs::Real theta = (std::atan2(y, x) >= 0) ? std::atan2(y, x) : std::atan2(y, x) + 2.0L * M_PI;

    return (2.0L / 3.0L) * std::pow(rho, -1.0L / 3.0L) * std::cos(theta / 3.0L);
}

/**
 * @brief Source.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real source(const pacs::Real &x, const pacs::Real &y) {
    return 0.0L;
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

/**
 * @brief Dirichlet boundary conditions, x derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real dirichlet_x(const pacs::Real &x, const pacs::Real &y) {
    return exact_x(x, y);
}

/**
 * @brief Dirichlet boundary conditions, y derivative.
 * 
 * @param x 
 * @param y 
 * @return pacs::Real 
 */
inline pacs::Real dirichlet_y(const pacs::Real &x, const pacs::Real &y) {
    return exact_y(x, y);
}