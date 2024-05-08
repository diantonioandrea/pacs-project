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

// Vector.
#include <Vector.hpp>

// OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// Containers.
#include <vector>

// Output.
#include <iostream>

// Assertions.
#include <cassert>

// Math.
#include <cmath>

// Copy.
#include <algorithm>

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
            Matrix(const std::size_t &rows, const std::size_t &columns, const std::vector<T> &elements): rows{rows}, columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert((rows > 0) && (columns > 0));
                assert(elements.size() == rows * columns);
                #endif

                this->elements.resize(rows * columns);
                std::ranges::copy(elements.begin(), elements.end(), this->elements.begin());
            }

            /**
             * @brief Copy constructor.
             * 
             * @param matrix 
             */
            Matrix(const Matrix &matrix): rows{matrix.rows}, columns{matrix.columns} {
                this->elements.resize(this->rows * this->columns);
                std::ranges::copy(matrix.elements.begin(), matrix.elements.end(), this->elements.begin());
            }
            
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

                this->elements.resize(matrix.elements.size());
                std::ranges::copy(matrix.elements.begin(), matrix.elements.end(), this->elements.begin());
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
                #ifndef NDEBUG // Integrity check.
                assert((j < this->rows) && (k < this->columns));
                #endif

                return this->elements[j * this->columns + k];
            }
            
            /**
             * @brief Call operator, returns a reference to the (i, j)-th element.
             * 
             * @param j 
             * @param k 
             * @return T& 
             */
            T &operator ()(const std::size_t &j, const std::size_t &k) {
                #ifndef NDEBUG // Integrity check.
                assert((j < this->rows) && (k < this->columns));
                #endif

                return this->elements[j * this->columns + k];
            }

            /**
             * @brief Returns the j-th row as a Vector.
             * 
             * @param j 
             * @return Vector<T> 
             */
            Vector<T> row(const std::size_t j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                #endif

                Vector<T> row{this->columns};

                for(std::size_t k = 0; k < this->columns; ++k)
                    row[k] = this->elements[j * this->columns + k];

                return row;
            }

            /**
             * @brief Sets the j-th row to the given scalar.
             * 
             * @param j 
             * @param scalar 
             */
            void row(const std::size_t j, const T &scalar) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                #endif

                for(std::size_t k = 0; k < this->columns; ++k)
                    this->elements[j * this->columns + k] = scalar;
            }

            /**
             * @brief Sets the j-th row to the given vector.
             * 
             * @param j 
             * @param vector 
             */
            void row(const std::size_t j, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                assert(vector.length == this->columns);
                #endif

                for(std::size_t k = 0; k < this->columns; ++k)
                    this->elements[j * this->columns + k] = vector[k];
            }

            /**
             * @brief Returns the k-th column as a Vector.
             * 
             * @param jk
             * @return Vector<T> 
             */
            Vector<T> column(const std::size_t &k) const {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->columns);
                #endif

                Vector<T> column{this->rows};

                for(std::size_t j = 0; j < this->rows; ++j)
                    column[j] = this->elements[j * this->columns + k];

                return column;
            }

            /**
             * @brief Sets the k-th column to the given scalar.
             * 
             * @param k 
             * @param scalar 
             */
            void column(const std::size_t &k, const T &scalar) {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->columns);
                #endif

                for(std::size_t j = 0; j < this->rows; ++j)
                    this->elements[j * this->columns + k] = scalar;
            }

            /**
             * @brief Sets the k-th column to the given vector.
             * 
             * @param k 
             * @param vector 
             */
            void column(const std::size_t &k, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->columns);
                assert(vector.length == this->rows);
                #endif

                for(std::size_t j = 0; j < this->rows; ++j)
                    this->elements[j * this->columns + k] = vector[j];
            }

            // SIZE.

            /**
             * @brief Returns the Matrix size.
             * 
             * @return std::size_t 
             */
            inline std::size_t size() const {
                return this->rows * this->columns;
            }

            // SHAPE.

            /**
             * @brief Returns the reshaped Matrix.
             * 
             * @param rows 
             * @param columns 
             * @return Matrix 
             */
            Matrix reshape(const std::size_t &rows, const std::size_t &columns) const {
                return Matrix{rows, columns, this->elements};
            }

            /**
             * @brief Returns the transpose matrix.
             * 
             * @return Matrix 
             */
            Matrix transpose() const {
                Matrix transpose{this->columns, this->rows};

                for(std::size_t j = 0; j < this->rows; ++j)
                    for(std::size_t k = 0; k < this->columns; ++k)
                        transpose.elements[k * this->rows + j] = this->elements[j * this->columns + k];

                return transpose;
            }

            // OPERATORS.

            /**
             * @brief Matrix scalar product.
             * 
             * @param scalar 
             * @return Matrix 
             */
            Matrix operator *(const T &scalar) const {
                Matrix result{*this};

                #pragma omp parallel for
                for(auto &element: result.elements)
                    element *= scalar;

                return result;
            }

            /**
             * @brief Friend Matrix scalar product.
             * 
             * @param scalar 
             * @param matrix 
             * @return Matrix 
             */
            friend Matrix operator *(const T &scalar, const Matrix &matrix) {
                Matrix result{matrix};

                #pragma omp parallel for
                for(auto &element: result.elements)
                    element *= scalar;

                return result;
            }

            /**
             * @brief Matrix scalar product and assignation.
             * 
             * @param scalar 
             * @return Matrix& 
             */
            Matrix &operator *=(const T scalar) {
                #pragma omp parallel for
                for(auto &element: this->elements)
                    element *= scalar;

                return *this;
            }

            /**
             * @brief Matrix sum.
             * 
             * @param matrix 
             * @return Matrix 
             */
            Matrix operator +(const Matrix &matrix) const {
                #ifndef NDEBUG // Integrity check.
                assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
                #endif

                Matrix result{*this};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->size(); ++j)
                    result.elements[j] += matrix.elements[j];

                return result;
            }

            /**
             * @brief Matrix sum and assignation.
             * 
             * @param matrix 
             * @return Matrix& 
             */
            Matrix &operator +=(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
                #endif

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->size(); ++j)
                    this->elements[j] += matrix.elements[j];

                return *this;
            }

            /**
             * @brief Matrix difference.
             * 
             * @param matrix 
             * @return Matrix 
             */
            Matrix operator -(const Matrix &matrix) const {
                #ifndef NDEBUG // Integrity check.
                assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
                #endif

                Matrix result{*this};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->size(); ++j)
                    result.elements[j] -= matrix.elements[j];

                return result;
            }

            /**
             * @brief Matrix difference and assignation.
             * 
             * @param matrix 
             * @return Matrix& 
             */
            Matrix &operator -=(const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
                #endif

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->size(); ++j)
                    this->elements[j] -= matrix.elements[j];

                return *this;
            }

            /**
             * @brief Matrix * Vector product.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> operator *(const Vector<T> &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->columns == vector.length);
                #endif

                Vector<T> result{this->rows};
                
                for(std::size_t j = 0; j < this->rows; ++j) {
                    T product = static_cast<T>(0);

                    #pragma omp parallel for reduction(+: product)
                    for(std::size_t k = 0; k < this->columns; ++k)
                        product += this->elements[j * this->columns + k] * vector[k];

                    result[j] = product;
                }

                return result;
            }

            /**
             * @brief Matrix * Matrix product.
             * 
             * @param matrix 
             * @return Matrix 
             */
            Matrix operator *(const Matrix &matrix) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->columns == matrix.rows);
                #endif

                Matrix result{this->rows, matrix.columns};
                
                #pragma omp parallel for collapse(2)
                for(std::size_t j = 0; j < this->rows; ++j)
                    for(std::size_t k = 0; k < matrix.columns; ++k) {
                        T product = static_cast<T>(0);

                        #pragma omp parallel for reduction (+: product)
                        for(std::size_t h = 0; h < this->columns; ++h)
                            product += this->elements[j * this->columns + h] * matrix.elements[h * matrix.columns + k];

                        result.elements[j * matrix.columns + k] = product;
                    }

                return result;
            }

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
                        ost << matrix.elements[j * matrix.columns + k] << " ";

                    if(j < matrix.rows - 1)
                        std::cout << std::endl;
                }

                return ost;
            }
    };

    // METHODS.

    /**
     * @brief Squashes a matrix to a vector.
     * 
     * @tparam T 
     * @param matrix 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> squash(const Matrix<T> &matrix) {
        Vector<T> result{matrix.rows * matrix.columns};

        #pragma omp parallel for collapse(2)
        for(std::size_t j = 0; j < matrix.rows; ++j)
            for(std::size_t k = 0; k < matrix.columns; ++k)
                result[j * matrix.columns + k] = matrix(j, k);

        return result;
    }

}

#endif