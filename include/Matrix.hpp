/**
 * @file Matrix.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MATRIX_PACS
#define MATRIX_PACS

// Type.
#include <Type.hpp>

// Containers.
#include <vector>
#include <array>
#include <map>

// Output.
#include <iostream>

// Assertions.
#include <cassert>

// Math.
#include <cmath>

// Zero tolerance.
#ifndef TOLERANCE_PACS
#define TOLERANCE_PACS 1E-10
#endif

namespace pacs {

    template<NumericType T>
    class Matrix {
        protected:

            // Elements.
            std::vector<T> elements;

        public:

            // Shape.
            const std::size_t rows;
            const std::size_t columns;

            // CONSTRUCTORS.

            // ...

            // READ AND WRITE.
            
            /**
             * @brief Const call operator, returns the (i, j)-th element.
             * 
             * @param j 
             * @param k 
             * @return T 
             */
            T operator ()(const std::size_t &j, const std::size_t &k) const {
                #ifndef NDEBUG // Out-of-bound check.
                assert((j < rows) && (k < columns));
                #endif

                return this->elements[j * this->rows + k];
            }
            
            /**
             * @brief Call operator, returns a reference to the (i, j)-th element.
             * 
             * @param j 
             * @param k 
             * @return T& 
             */
            T &operator ()(const std::size_t &j, const std::size_t &k) {
                #ifndef NDEBUG // Out-of-bound check.
                assert((j < rows) && (k < columns));
                #endif

                return this->elements[j * this->rows + k];
            }

            // SHAPE

            Matrix reshape(const std::size_t &rows, const std::size_t &columns) const {
                return Matrix{rows, colums, this->elements};
            }
    };

}

#endif