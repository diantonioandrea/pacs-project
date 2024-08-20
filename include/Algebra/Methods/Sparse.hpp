/**
 * @file Sparse.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SPARSE_METHODS_PACS
#define SPARSE_METHODS_PACS

#include "../Sparse.hpp"

namespace pacs {

    /**
     * @brief Product of diagonal elements.
     * 
     * @tparam T Numeric type.
     * @param sparse Sparse matrix.
     * @return T
     */
    template<NumericType T>
    T mtrace(const Sparse<T> &sparse) {
        #ifndef NDEBUG // Integrity check.
        assert(sparse.rows == sparse.columns);
        #endif

        T product = static_cast<T>(1);

        for(std::size_t j = 0; j < sparse.rows; ++j)
            product *= sparse(j, j);

        return product;
    }
    
}

#endif