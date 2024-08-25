/**
 * @file Fit.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <cassert>

namespace pacs {

    /**
     * @brief Polynomial fit.
     * 
     * @param x X.
     * @param y Y.
     * @param p Algorithm order.
     * @return Vector<Real> 
     */
    Vector<Real> polyfit(const Vector<Real> &x, const Vector<Real> &y, const std::size_t &p) {
        #ifndef NDEBUG // Integrity check.
        assert(p > 0);
        assert(x.length == y.length);
        #endif

        // X.
        Matrix<Real> X{x.length, p + 1};

        // Building X.
        for(std::size_t j = 0; j < p + 1; ++j) {
            X.column(j, Vector<Real>{x.length, 1.0});

            for(std::size_t k = 0; k < j; ++k) {
                X.column(j, X.column(j) * x);
            }
        }

        // Solution.
        return solve(X.transpose() * X, X.transpose() * y, QRD);
    }

}