/**
 * @file Legendre.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LEGENDRE_PACS
#define LEGENDRE_PACS

// Vectors.
#include <Vector.hpp>

namespace pacs {

    // Legendre polynomials.

    Vector<double> legendre(const Vector<double> &, const std::size_t &);
    Vector<double> grad_legendre(const Vector<double> &, const std::size_t &);

}

#endif