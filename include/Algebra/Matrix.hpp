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

#include "Vector.hpp"

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>

namespace pacs {

    /**
     * @brief Matrix structure.
     * 
     * @tparam T 
     */
    template<NumericType T>
    struct Matrix {

        // Shape.
        const std::size_t rows;
        const std::size_t columns;

        // Elements.
        std::vector<T> elements;

        // CONSTRUCTORS.

        /**
         * @brief Constructs a new empty Matrix.
         * 
         * @param rows 
         * @param columns 
         */
        inline Matrix(const std::size_t &rows, const std::size_t &columns): rows{rows}, columns{columns}, elements(rows * columns, static_cast<T>(0)) {
            #ifndef NDEBUG // Integrity check.
            assert((rows > 0) && (columns > 0));
            #endif
        }

        /**
         * @brief Constructs a new Matrix from a given std::vector.
         * 
         * @param rows 
         * @param columns 
         * @param elements 
         */
        inline Matrix(const std::size_t &rows, const std::size_t &columns, const std::vector<T> &elements): rows{rows}, columns{columns}, elements(elements.begin(), elements.end()) {
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
        inline Matrix(const Matrix &matrix): rows{matrix.rows}, columns{matrix.columns}, elements(matrix.elements.begin(), matrix.elements.end()) {}
        
        /**
         * @brief Copy operator.
         * 
         * @param matrix 
         * @return Matrix& 
         */
        inline Matrix &operator =(const Matrix &matrix) {
            #ifndef NDEBUG
            assert((this->rows == matrix.rows) && (this->columns == matrix.columns));
            #endif

            #ifdef PARALLEL
            std::copy(POLICY, matrix.elements.begin(), matrix.elements.end(), this->elements.begin());
            #else
            std::copy(matrix.elements.begin(), matrix.elements.end(), this->elements.begin());
            #endif

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
        inline T operator ()(const std::size_t &j, const std::size_t &k) const {
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
        inline T &operator ()(const std::size_t &j, const std::size_t &k) {
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

            #ifdef PARALLEL
            std::copy(POLICY, this->elements.begin() + j * this->columns, this->elements.begin() + (j + 1) * this->columns, row.elements.begin());
            #else
            std::copy(this->elements.begin() + j * this->columns, this->elements.begin() + (j + 1) * this->columns, row.elements.begin());
            #endif

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

            #ifdef PARALLEL
            std::for_each_n(POLICY, this->elements.begin() + j * this->columns, this->columns, [scalar](auto &element){ element = scalar; });
            #else
            std::for_each_n(this->elements.begin() + j * this->columns, this->columns, [scalar](auto &element){ element = scalar; });
            #endif
        }

        /**
         * @brief Sets the j-th row to the given Vector.
         * 
         * @param j 
         * @param vector 
         */
        void row(const std::size_t j, const Vector<T> &vector) {
            #ifndef NDEBUG // Integrity check.
            assert(j < this->rows);
            assert(vector.length == this->columns);
            #endif

            #ifdef PARALLEL
            std::copy(POLICY, vector.elements.begin(), vector.elements.end(), this->elements.begin() + j * this->columns);
            #else
            std::copy(vector.elements.begin(), vector.elements.end(), this->elements.begin() + j * this->columns);
            #endif
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
                this->elements[j * this->columns + k] = vector.elements[j];
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
        inline Matrix reshape(const std::size_t &rows, const std::size_t &columns) const {
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
        
        /**
         * @brief Returns the diagonal.
         * 
         * @return Matrix 
         */
        Matrix diagonal() const {
            #ifndef NDEBUG // Integrity check.
            assert(this->rows == this->columns);
            #endif

            Matrix diagonal{this->rows, this->columns};

            for(std::size_t j = 0; j < this->rows; ++j)
                diagonal.elements[j * (this->columns + 1)] = this->elements[j * (this->columns + 1)];

            return diagonal;
        }

        /**
         * @brief Returns the lower triangular part of the Matrix.
         * 
         * @return Matrix 
         */
        Matrix lower() const {
            #ifndef NDEBUG // Integrity check.
            assert(this->rows == this->columns);
            #endif

            Matrix lower{this->rows, this->columns};

            for(std::size_t j = 0; j < this->rows; ++j)
                #ifdef PARALLEL
                std::copy(POLICY, this->elements.begin() + j * this->columns, this->elements.begin() + j * this->columns + j, lower.elements.begin() + j * this->columns);
                #else
                std::copy(this->elements.begin() + j * this->columns, this->elements.begin() + j * this->columns + j, lower.elements.begin() + j * this->columns);
                #endif

            return lower;
        }

        /**
         * @brief Returns the upper triangular part of the Matrix.
         * 
         * @return Matrix 
         */
        Matrix upper() const {
            #ifndef NDEBUG // Integrity check.
            assert(this->rows == this->columns);
            #endif

            Matrix upper{this->rows, this->columns};

            for(std::size_t j = 0; j < this->rows; ++j)
                #ifdef PARALLEL
                std::copy(POLICY, this->elements.begin() + j * this->columns + j + 1, this->elements.begin() + (j + 1) * this->columns, upper.elements.begin() + j * this->columns + j + 1);
                #else
                std::copy(this->elements.begin() + j * this->columns + j + 1, this->elements.begin() + (j + 1) * this->columns, upper.elements.begin() + j * this->columns + j + 1);
                #endif

            return upper;
        }

        // OPERATIONS.

        /**
         * @brief Matrix unary +.
         * 
         * @return Matrix 
         */
        Matrix operator +() const {
            return *this;
        }

        /**
         * @brief Matrix unary -.
         * 
         * @return Matrix 
         */
        Matrix operator -() const {
            Matrix result{this->rows, this->columns};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [](const auto &element){ return -element; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [](const auto &element){ return -element; });
            #endif

            return result;
        }

        /**
         * @brief Matrix scalar product.
         * 
         * @param scalar 
         * @return Matrix 
         */
        Matrix operator *(const T &scalar) const {
            Matrix result{this->rows, this->columns};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

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
            Matrix result{matrix.rows, matrix.columns};

            #ifdef PARALLEL
            std::transform(POLICY, matrix.elements.begin(), matrix.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(matrix.elements.begin(), matrix.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

            return result;
        }

        /**
         * @brief Matrix scalar product and assignation.
         * 
         * @param scalar 
         * @return Matrix& 
         */
        Matrix &operator *=(const T scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

            return *this;
        }

        /**
         * @brief Matrix scalar division.
         * 
         * @param scalar 
         * @return Matrix 
         */
        Matrix operator /(const T &scalar) const {
            Matrix result{this->rows, this->columns};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #endif

            return result;
        }

        /**
         * @brief Matrix scalar division and assignation.
         * 
         * @param scalar 
         * @return Matrix& 
         */
        Matrix &operator /=(const T scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #endif

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

            Matrix result{this->rows, this->columns};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), matrix.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), matrix.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #endif

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

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), matrix.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), matrix.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #endif

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

            Matrix result{this->rows, this->columns};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), matrix.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), matrix.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #endif

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

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), matrix.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), matrix.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #endif

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

            for(std::size_t j = 0; j < this->rows; ++j)
                #ifdef PARALLEL
                result.elements[j] = std::transform_reduce(POLICY, vector.elements.begin(), vector.elements.end(), this->elements.begin() + j * this->columns, static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * second; });
                #else 
                result.elements[j] = std::inner_product(vector.elements.begin(), vector.elements.end(), this->elements.begin() + j * this->columns, static_cast<T>(0));
                #endif

            return result;
        }

        /**
         * @brief Friend Vector * Matrix product.
         * 
         * @param vector 
         * @param matrix 
         * @return Vector<T> 
         */
        friend Vector<T> operator *(const Vector<T> &vector, const Matrix &matrix) {
            #ifndef NDEBUG // Integrity check.
            assert(vector.length == matrix.rows);
            #endif

            Vector<T> result{matrix.columns};

            for(std::size_t j = 0; j < matrix.columns; ++j)
                for(std::size_t k = 0; k < matrix.rows; ++k)
                    result.elements[j] += vector.elements[k] * matrix.elements[k * matrix.rows + j];

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
            
            for(std::size_t j = 0; j < this->rows; ++j)
                for(std::size_t k = 0; k < matrix.columns; ++k)
                    for(std::size_t h = 0; h < this->columns; ++h)
                        result.elements[j * result.columns + k] += this->elements[j * this->columns + h] * matrix.elements[h * matrix.columns + k];

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

}

#endif