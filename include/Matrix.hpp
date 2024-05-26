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
            Matrix(const std::size_t &rows, const std::size_t &columns): elements(rows * columns, static_cast<T>(0)), rows{rows}, columns{columns} {
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
                    for(std::size_t k = 0; k < j; ++k)
                        lower.elements[j * this->columns + k] = this->elements[j * this->columns + k];

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
                    for(std::size_t k = this->columns - 1; k > j; --k)
                        upper.elements[j * this->columns + k] = this->elements[j * this->columns + k];

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
                Matrix result{*this};

                for(auto &element: result.elements)
                    element = -element;

                return result;
            }

            /**
             * @brief Matrix scalar product.
             * 
             * @param scalar 
             * @return Matrix 
             */
            Matrix operator *(const T &scalar) const {
                Matrix result{*this};

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
                for(auto &element: this->elements)
                    element *= scalar;

                return *this;
            }

            /**
             * @brief Matrix scalar division.
             * 
             * @param scalar 
             * @return Matrix 
             */
            Matrix operator /(const T &scalar) const {
                Matrix result{*this};

                for(auto &element: result.elements)
                    element /= scalar;

                return result;
            }

            /**
             * @brief Matrix scalar division and assignation.
             * 
             * @param scalar 
             * @return Matrix& 
             */
            Matrix &operator /=(const T scalar) {
                for(auto &element: this->elements)
                    element /= scalar;

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

                    for(std::size_t k = 0; k < this->columns; ++k)
                        product += this->elements[j * this->columns + k] * vector[k];

                    result[j] = product;
                }

                return result;
            }

            friend Vector<T> operator *(const Vector<T> &vector, const Matrix &matrix) {
                #ifndef NDEBUG // Integrity check.
                assert(vector.length == matrix.rows);
                #endif

                Vector<T> result{matrix.columns};

                for(std::size_t j = 0; j < matrix.columns; ++j) {
                    T product = static_cast<T>(0);

                    for(std::size_t k = 0; k < matrix.rows; ++k)
                        product += matrix.elements[k * matrix.columns + j];

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
                
                for(std::size_t j = 0; j < this->rows; ++j)
                    for(std::size_t k = 0; k < matrix.columns; ++k)
                        for(std::size_t h = 0; h < this->columns; ++h)
                            result.elements[j * result.columns + k] += this->elements[j * this->columns + h] * matrix.elements[h * matrix.columns + k];

                return result;
            }

            // DECOMPOSITIONS.

            /**
             * @brief Returns the LU decomposition of the Matrix.
             * 
             * @return std::array<Matrix, 2> 
             */
            std::array<Matrix, 2> LU() const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                #endif

                // Target.
                Matrix target{*this};

                // LU decomposition
                Matrix L{this->rows, this->columns};
                Matrix U{this->rows, this->columns};

                for (std::size_t j = 0; j < this->columns; ++j) {
                    L(j, j) = static_cast<T>(1);

                    // Compute elements of U
                    for (std::size_t i = 0; i <= j; ++i) {
                        T sum = static_cast<T>(0);

                        for (std::size_t k = 0; k < i; ++k)
                            sum += L(i, k) * U(k, j);

                        U(i, j) = target(i, j) - sum;
                    }

                    // Compute elements of L
                    for (std::size_t i = j + 1; i < this->rows; ++i) {
                        T sum = static_cast<T>(0);

                        for (std::size_t k = 0; k < j; ++k)
                            sum += L(i, k) * U(k, j);

                        L(i, j) = (target(i, j) - sum) / U(j, j);
                    }
                }

                return {L, U};
            }

            /**
             * @brief Returns the QR decomposition of the Matrix.
             * 
             * @return std::array<Matrix, 2> 
             */
            std::array<Matrix, 2> QR() const {

                // Identity matrix.
                Matrix I{this->rows, this->rows};

                for(std::size_t j = 0; j < this->rows; ++j)
                    I.elements[j * (this->rows + 1)] = static_cast<T>(1);

                // QR decomposition.
                Matrix Q{I};
                Matrix R{*this};

                // Algorithm, real case.
                for(std::size_t j = 0; j < ((this->rows > this->columns) ? this->columns : this->rows); ++j) {

                    // Householder vector.
                    Vector<T> vector{this->rows - j};

                    for(std::size_t k = 0; k < this->rows - j; ++k)
                        vector[k] = R.elements[(j + k) * this->columns + j];

                    vector[0] += (vector[0] > 0 ? static_cast<T>(1) : static_cast<T>(-1)) * norm(vector);
                    vector /= norm(vector);

                    // Householder matrix.
                    Matrix H{I};

                    for(std::size_t k = 0; k < this->rows - j; ++k)
                        for(std::size_t l = 0; l < this->rows - j; ++l)
                            H.elements[(j + k) * this->rows + (j + l)] -= 2 * vector[k] * vector[l];

                    // QR.
                    R = H * R;
                    Q = Q * H.transpose();
                }

                return {Q, R};
            }

            // LINEAR.

            /**
             * @brief Solves a linear system in the form of Ax = b.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> solve(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == vector.length);
                #endif

                if(this->rows == this->columns)
                    return lu_solver(vector);

                return qr_solver(vector);
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

        private:

            // SOLVERS.

            /**
             * @brief LU solver.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> lu_solver(const Vector<T> &vector) const {

                // LU decomposition.
                auto [L, U] = this->LU();

                // Solves Ly = b using forward substitution.
                Vector<T> y{this->rows};

                for (std::size_t i = 0; i < this->rows; ++i) {
                    T sum = static_cast<T>(0);

                    for (std::size_t j = 0; j < i; ++j)
                        sum += L(i, j) * y[j];

                    y[i] = (vector[i] - sum) / L(i, i);
                }

                // Solves Ux = y using backward substitution.
                Vector<T> x{this->columns};

                for (std::size_t i = this->columns; i > 0; --i) {
                    T sum = static_cast<T>(0);

                    for (std::size_t j = i; j < this->columns; ++j)
                        sum += U(i - 1, j) * x[j];

                    x[i - 1] = (y[i - 1] - sum) / U(i - 1, i - 1);
                }

                return x;
            }
            
            /**
             * @brief QR solver.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> qr_solver(const Vector<T> &vector) const {
                
                // QR decomposition.
                auto [Q, R] = this->QR();

                // Solves Rx = QTb using backward substitution.
                Vector<T> x{this->columns};
                Vector<T> b{Q.transpose() * vector};

                for (std::size_t i = this->columns; i > 0; --i) {
                    T sum = static_cast<T>(0);

                    for (std::size_t j = i; j < this->columns; ++j)
                        sum += R(i - 1, j) * x[j];

                    x[i - 1] = (b[i - 1] - sum) / R(i - 1, i - 1);
                }

                return x;
            }
    };

    // METHODS.

    template<NumericType T>
    inline Matrix<T> identity(const std::size_t &size) {
        #ifndef NDEBUG // Integrity check.
        assert(size > 0);
        #endif

        Matrix<T> I{size, size};

        for(std::size_t j = 0; j < size; ++j)
            I(j, j) = static_cast<T>(1);

        return I;
    }

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

        for(std::size_t j = 0; j < matrix.rows; ++j)
            for(std::size_t k = 0; k < matrix.columns; ++k)
                result[j * matrix.columns + k] = matrix(j, k);

        return result;
    }

    /**
     * @brief Multiplicative trace.
     * 
     * @tparam T 
     * @param matrix 
     * @return T 
     */
    template<NumericType T>
    T mtrace(const Matrix<T> &matrix) {
        #ifndef NDEBUG // Integrity check.
        assert(matrix.rows == matrix.columns);
        #endif

        T product = static_cast<T>(1);

        for(std::size_t j = 0; j < matrix.rows; ++j)
            product *= matrix(j, j);

        return product;
    }

}

#endif