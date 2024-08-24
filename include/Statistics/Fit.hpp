/**
 * @file Fit.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FIT_PACS
#define FIT_PACS

#include "../Base.hpp"
#include "../Algebra.hpp"

namespace pacs {

    // Polynomial fit.

    Vector<Real> polyfit(const Vector<Real> &, const Vector<Real> &, const std::size_t &);

}

#endif