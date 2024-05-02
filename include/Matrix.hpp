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

// Output.
#include <iostream>

// Assertions.
#include <cassert>

// Math.
#include <cmath>

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

            /**
             * @brief Constructs a new empty Matrix.
             * 
             * @param rows 
             * @param columns 
             */
            Matrix(const std::size_t &rows, const std::size_t &columns): rows{rows}, columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert((rows > 0) && (columns > 0));
                #endif

                this->elements.resize(rows * columns, static_cast<T>(0));
            }

            /**
             * @brief Constructs a new Matrix from a given std::vector.
             * 
             * @param rows 
             * @param columns 
             * @param elements 
             */
            Matrix(const std::size_t &rows, const std::size_t &columns, const std::vector<T> &elements): rows{rows}, columns{columns}, elements{elements} {
                #ifndef NDEBUG // Integrity check.
                assert((rows > 0) && (columns > 0));
                assert(elements.size() == rows * columns);
                #endif
            }

            /**
             * @brief Copy constructor.
             * 
             * @param matrix 
             */
            Matrix(const Matrix &matrix): rows{matrix.rows}, columns{matrix.columns}, elements{matrix.elements} {}
            
            /**
             * @brief Copy operator.
             * 
             * @param matrix 
             * @return Matrix& 
             */
            Matrix &operator =(const Matrix &matrix) {
                #ifndef NDEBUG
                assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
                #endif

                this->elements = matrix.elements;
                return *this;
            }

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

            // SHAPE.

            Matrix reshape(const std::size_t &rows, const std::size_t &columns) const {
                return Matrix{rows, columns, this->elements};
            }

            // OPERATORS.

            // ...

            // OUTPUT.

            /**
             * @brief Matrix output.
             * 
             * @param ost 
             * @param matrix 
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Matrix &matrix) {
                for(std::size_t j = 0; j < matrix.rows; ++j) {
                    for(std::size_t k = 0; k < matrix.columns; ++k)
                        ost << matrix.elements[j * matrix.rows + k] << " ";

                    if(j < matrix.columns - 1)
                        std::cout << std::endl;
                }

                return ost;
            }
    };

}

#endif