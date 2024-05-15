/**
 * @file Sparse.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SPARSE_PACS
#define SPARSE_PACS

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
#include <array>
#include <map>

// Output.
#include <iostream>

// Assertions.
#include <cassert>

// Math.
#include <cmath>

// Matrices.
#include <Matrix.hpp>

// Copy.
#include <algorithm>

// Algebra tolerance.
#ifndef ALGEBRA_TOLERANCE
#define ALGEBRA_TOLERANCE 1E-4
#endif

// Algebra iterations limit.
#ifndef ALGEBRA_ITER_MAX
#define ALGEBRA_ITER_MAX 1E5
#endif

namespace pacs {

    // Solvers.
    enum Solver {Conjugate, Descent, Gauss, Kaczmarz, RandomKaczmarz};
    
    /**
     * @brief Sparse matrix class.
     * 
     * @tparam T Matrix' type.
     */
    template<NumericType T>
    class Sparse {
        protected:

            // Compression flag.
            bool compressed = false;

            // COOmap dynamic storage format.
            mutable std::map<std::array<std::size_t, 2>, T> elements;

            // CSR compressed storage format.
            std::vector<std::size_t> inner;
            std::vector<std::size_t> outer;
            std::vector<T> values;

        public:

            // Shape.
            const std::size_t rows;
            const std::size_t columns;

            // CONSTRUCTORS.

            /**
             * @brief Constructs a new empty Sparse matrix.
             * 
             * @param rows 
             * @param columns 
             */
            Sparse(const std::size_t &rows, const std::size_t &columns): rows{rows}, columns{columns} {
                #ifndef NDEBUG // Integrity check.
                assert((rows > 0) && (columns > 0));
                #endif
            }

            /**
             * @brief Constructs a new Sparse matrix from a given std::map.
             * 
             * @param rows 
             * @param columns 
             * @param elements 
             */
            Sparse(const std::size_t &rows, const std::size_t &columns, const std::map<std::array<std::size_t, 2>, T> &elements): rows{rows}, columns{columns}, elements{elements} {
                #ifndef NDEBUG // Integrity checks.
                assert((rows > 0) && (columns > 0));

                for(const auto &[key, value]: elements)
                    assert((key[0] < rows) && (key[1] < columns));

                #endif
            }

            /**
             * @brief Constructs a new Sparse matrix from given inner, outer and values vectors.
             *
             * @param rows
             * @param columns
             * @param elements
             */
            Sparse(const std::size_t &rows, const std::size_t &columns, const std::vector<std::size_t> &inner, const std::vector<std::size_t> &outer, const std::vector<T> &values):
            rows{rows}, columns{columns}, compressed{true}, inner{inner}, outer{outer}, values{values} {
                #ifndef NDEBUG // Integrity checks.
                assert((rows > 0) && (columns > 0));
                assert(inner.size() == rows + 1);
                assert(outer.size() = values.size());

                for(std::size_t j = 1; j < inner.size(); ++j) {
                    assert(inner[j - 1] < values.size());
                    assert(inner[j] < values.size());
                    assert(inner[j - 1] <= inner[j]);
                }

                for(std::size_t j = 1; j < outer.size(); ++j) {
                    assert(outer[j - 1] < columns);
                    assert(outer[j] < columns);
                    assert(outer[j - 1] < outer[j]);
                }
                #endif
            }

            /**
             * @brief Copy constructor.
             *
             * @param sparse
             */
            Sparse(const Sparse &sparse): compressed{sparse.compressed}, rows{sparse.rows}, columns{sparse.columns} {
                if(!(sparse.compressed))
                    this->elements = sparse.elements;
                else {
                    this->inner.resize(sparse.inner.size());
                    this->outer.resize(sparse.outer.size());
                    this->values.resize(sparse.values.size());

                    std::ranges::copy(sparse.inner.begin(), sparse.inner.end(), this->inner.begin());
                    std::ranges::copy(sparse.outer.begin(), sparse.outer.end(), this->outer.begin());
                    std::ranges::copy(sparse.values.begin(), sparse.values.end(), this->values.begin());
                }
            }

            /**
             * @brief Copy operator.
             *
             * @param sparse
             * @return Sparse&
             */
            Sparse &operator =(const Sparse &sparse) {
                #ifndef NDEBUG
                assert((this->rows == sparse.rows) && (this->columns == sparse.columns));
                #endif

                this->compressed = sparse.compressed;

                this->elements.clear();
                this->inner.clear();
                this->outer.clear();
                this->values.clear();

                if(!(sparse.compressed))
                    this->elements = sparse.elements;
                else {
                    this->inner.resize(sparse.inner.size());
                    this->outer.resize(sparse.outer.size());
                    this->values.resize(sparse.values.size());

                    std::ranges::copy(sparse.inner.begin(), sparse.inner.end(), this->inner.begin());
                    std::ranges::copy(sparse.outer.begin(), sparse.outer.end(), this->outer.begin());
                    std::ranges::copy(sparse.values.begin(), sparse.values.end(), this->values.begin());
                }

                return *this;
            }

            // READ.

            /**
             * @brief Const call operator, returns the (i, j)-th element if present.
             *
             * @param j
             * @param k
             * @return T
             */
            T operator ()(const std::size_t &j, const std::size_t &k) const {
                #ifndef NDEBUG // Out-of-bound check.
                assert((j < this->rows) && (k < this->columns));
                #endif

                // Checks for the value inside elements, otherwise returns static_cast<T>(0).
                if(!(this->compressed))
                    return this->elements.contains({j, k}) ? this->elements[{j, k}] : static_cast<T>(0);

                // Looks for the value on compressed Matrix.
                for(std::size_t i = this->inner[j]; i < this->inner[j + 1]; ++i)
                    if(k == this->outer[i])
                        return this->values[i];

                // Default return.
                return static_cast<T>(0);
            }

            /**
             * @brief Const call operator, returns a sub-matrix.
             * 
             * @param J 
             * @param K 
             * @return Matrix<T> 
             */
            Matrix<T> operator ()(const std::vector<std::size_t> &J, const std::vector<std::size_t> &K) const {
                #ifndef NDEBUG // Out-of-bound check.
                for(std::size_t j = 0; j < J.size(); ++j)
                    for(std::size_t k = 0; k < K.size(); ++k)
                        assert((J[j] < this->rows) && (K[k] < this->columns));
                #endif

                Matrix<T> matrix{J.size(), K.size()};

                for(std::size_t j = 0; j < J.size(); ++j)
                    for(std::size_t k = 0; k < K.size(); ++k)
                        matrix(j, k) = (*this)(J[j], K[k]);
                
                return matrix;
            }

            // INSERT.

            /**
             * @brief Inserts a new element.
             * 
             * @param j 
             * @param k 
             * @param element 
             */
            void insert(const std::size_t &j, const std::size_t &k, const T &element) {
                #ifndef NDEBUG // Integrity check.
                assert((j < this->rows) && (k < this->columns));
                assert(!(this->compressed));
                
                if(std::abs(element) > TOLERANCE)
                    this->elements[{j, k}] = element;
                #else
                this->elements[{j, k}] = element;
                #endif
            }

            /**
             * @brief Inserts a new matrix of elements.
             * 
             * @param j 
             * @param k 
             * @param elements 
             */
            void insert(const std::vector<std::size_t> &J, const std::vector<std::size_t> &K, const Matrix<T> &elements) {
                #ifndef NDEBUG // Integrity checks.
                for(std::size_t j = 0; j < J.size(); ++j)
                    assert((J[j] < this->rows) && (j < elements.rows));
                for(std::size_t k = 0; k < K.size(); ++k)
                    assert((K[k] < this->columns) && (k < elements.columns));
                #endif

                for(std::size_t j = 0; j < J.size(); ++j)
                    for(std::size_t k = 0; k < K.size(); ++k)
                        if(std::abs(elements(j, k)) > TOLERANCE)
                            this->elements[{J[j], K[k]}] = elements(j, k);
            }

            // ADD.

            /**
             * @brief Adds a new element. Slower than an insert for simple creation.
             * 
             * @param j 
             * @param k 
             * @param element 
             */
            void add(const std::size_t &j, const std::size_t &k, const T &element) {
                #ifndef NDEBUG // Integrity check.
                assert((j < this->rows) && (k < this->columns));
                assert(!(this->compressed));
                
                if(std::abs(element) > TOLERANCE) {
                    if(this->elements.contains({j, k}))
                        this->elements[{j, k}] += element;
                    else
                        this->elements[{j, k}] = element;
                }
                #else
                this->elements[{j, k}] += element;
                #endif
            }

            /**
             * @brief Adds a new matrix of elements.
             * 
             * @param j 
             * @param k 
             * @param elements 
             */
            void add(const std::vector<std::size_t> &J, const std::vector<std::size_t> &K, const Matrix<T> &elements) {
                #ifndef NDEBUG // Integrity checks.
                for(std::size_t j = 0; j < J.size(); ++j)
                    assert((J[j] < this->rows) && (j < elements.rows));
                for(std::size_t k = 0; k < K.size(); ++k)
                    assert((K[k] < this->columns) && (k < elements.columns));
                #endif

                for(std::size_t j = 0; j < J.size(); ++j)
                    for(std::size_t k = 0; k < K.size(); ++k)
                        if(std::abs(elements(j, k)) > TOLERANCE) {
                            if(this->elements.contains({J[j], K[k]}))
                                this->elements[{J[j], K[k]}] += elements(j, k);
                            else
                                this->elements[{J[j], K[k]}] = elements(j, k);
                        }
            }

            // ROW AND COLUMN.

            /**
             * @brief Returns the j-th row as a Vector.
             * 
             * @param j 
             * @return Vector<T> 
             */
            Vector<T> row(const std::size_t &j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                #endif

                Vector<T> row{this->columns};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements) {
                        if(key[0] > j)
                            break;

                        if(key[0] == j)
                            row[key[1]] = element;
                    }
                else
                    for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                        row[this->outer[k]] = this->values[k];

                return row;
            }

            /**
             * @brief Sets the j-th row to the given scalar.
             * 
             * @param j 
             * @param scalar 
             */
            void row(const std::size_t &j, const T &scalar) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                assert(!this->compressed);
                #endif

                if(std::abs(scalar) <= TOLERANCE)
                    return;

                for(std::size_t k = 0; k < this->columns; ++k)
                    this->elements[{j, k}] = scalar;
            }

            /**
             * @brief Sets the j-th row to the given Vector.
             * 
             * @param j 
             * @param vector 
             */
            void row(const std::size_t &j, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->rows);
                assert(vector.length == this->columns);
                assert(!this->compressed);
                #endif

                for(std::size_t k = 0; k < this->columns; ++k)
                    this->insert(j, k, vector[k]);
            }

            /**
             * @brief Return the k-th column as a Vector.
             * 
             * @param k 
             * @return Vector<T> 
             */
            Vector<T> column(const std::size_t &k) const {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->columns);
                #endif

                Vector<T> column{this->rows};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements) {
                        if(key[1] == k)
                            column[key[0]] = element;
                    }
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t h = this->inner[j]; h < this->inner[j + 1]; ++h)
                            if(this->outer[h] == k) {
                                column[j] = this->values[h];
                                break;
                            }

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
                assert(!this->compressed);
                #endif

                if(std::abs(scalar) <= TOLERANCE)
                    return;

                for(std::size_t j = 0; j < this->rows; ++j)
                    this->elements[{j, k}] = scalar;
            }

            /**
             * @brief Sets the k-th column to the given Vector.
             * 
             * @param k 
             * @param vector 
             */
            void column(const std::size_t &k, const Vector<T> &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(k < this->columns);
                assert(vector.length == this->rows);
                assert(!this->compressed);
                #endif

                for(std::size_t j = 0; j < this->rows; ++j)
                    this->insert(j, k, vector[j]);
            }

            // SHAPE.

            /**
             * @brief Returns the reshaped Sparse matrix.
             * 
             * @param rows 
             * @param columns 
             * @return Sparse 
             */
            Sparse reshape(const std::size_t &rows, const std::size_t &columns) const {
                if(!(this->compressed))
                    return Sparse{rows, columns, this->elements};

                return Sparse{rows, columns, this->inner, this->outer, this->values};
            }

            /**
             * @brief Returns the transpose Sparse matrix.
             * 
             * @return Sparse 
             */
            Sparse transpose() const {
                Sparse transpose{this->columns, this->rows};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements)
                        transpose.elements[{key[1], key[0]}] = element;
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            transpose.elements[{this->outer[k], j}] = this->values[k];

                return transpose;
            }

            /**
             * @brief Returns the diagonal.
             * 
             * @return Sparse 
             */
            Sparse diagonal() const {
                #ifndef NDEBUG // Integrity check.
                assert(this->rows == this->columns);
                #endif

                Sparse diagonal{this->rows, this->columns};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements) {
                        if(key[0] == key[1])
                            diagonal.elements[key] = element;
                    }
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j == this->outer[k]) {
                                diagonal.elements[{j, this->outer[k]}] = this->values[k];
                                break;
                            }
                
                return diagonal;
            }

            /**
             * @brief Returns the lower triangular part of the Sparse matrix.
             * 
             * @return Sparse 
             */
            Sparse lower() const {
                #ifndef NDEBUG // Integrity check.
                assert(this->rows == this->columns);
                #endif

                Sparse lower{this->rows, this->columns};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements) {
                        if(key[0] > key[1])
                            lower.elements[key] = element;
                    }
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j > this->outer[k])
                                lower.elements[{j, this->outer[k]}] = this->values[k];
                
                return lower;
            }

            /**
             * @brief Returns the upper triangular part of the Sparse matrix.
             * 
             * @return Sparse 
             */
            Sparse upper() const {
                #ifndef NDEBUG // Integrity check.
                assert(this->rows == this->columns);
                #endif

                Sparse upper{this->rows, this->columns};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements) {
                        if(key[0] < key[1])
                            upper.elements[key] = element;
                    }
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j < this->outer[k])
                                upper.elements[{j, this->outer[k]}] = this->values[k];
                
                return upper;
            }

            // SHAPE CHECKS.

            /**
             * @brief Checks whether the Sparse matrix is diagonal.
             * 
             * @return true 
             * @return false 
             */
            bool is_diagonal() const {
                if(!(this->compressed)) {
                    for(const auto &[key, element]: this->elements)
                        if(key[0] != key[1])
                            return false;
                } else {
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j != this->outer[k])
                                return false;
                }

                return true;
            }

            /**
             * @brief Checks whether the Sparse matrix is lower triangular.
             * 
             * @return true 
             * @return false 
             */
            bool is_lower() const {
                if(!(this->compressed)) {
                    for(const auto &[key, element]: this->elements)
                        if(key[0] < key[1])
                            return false;
                } else {
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j < this->outer[k])
                                return false;
                }

                return true;
            }

            /**
             * @brief Checks whether the Sparse matrix is upper triangular.
             * 
             * @return true 
             * @return false 
             */
            bool is_upper() const {
                if(!(this->compressed)) {
                    for(const auto &[key, element]: this->elements)
                        if(key[0] > key[1])
                            return false;
                } else {
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            if(j > this->outer[k])
                                return false;
                }

                return true;
            }

            // COMPRESSION.

            /**
             * @brief Compresses an uncompressed Sparse matrix.
             *
             */
            void compress() {
                if(this->compressed)
                    return;

                std::size_t index = 0;
                std::array<std::size_t, 2> current{0, 0};
                std::array<std::size_t, 2> next{1, 0};

                this->inner.resize(this->rows + 1);
                this->inner[0] = index;

                // Compression.
                for(std::size_t j = 1; j < this->rows + 1; ++j) {
                    for(auto it = this->elements.lower_bound(current); (*it).first < (*(this->elements.lower_bound(next))).first; ++it) {
                        auto [key, value] = (*it);

                        if(std::abs(value) > TOLERANCE) {
                            this->outer.emplace_back(key[1]);
                            this->values.emplace_back(value);
                            ++index;
                        }
                    }

                    this->inner[j] = index;
                    ++current[0];
                    ++next[0];
                }

                this->compressed = true;
                this->elements.clear();
            }

            /**
             * @brief Uncompresses a compressed Sparse matrix.
             *
             */
            void uncompress() {
                if(!(this->compressed))
                    return;

                // Uncompression.
                for(std::size_t j = 0; j < this->inner.size() - 1; ++j)
                    for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                        this->elements[{j, this->outer[k]}] = this->values[k];

                this->compressed = false;
                this->inner.clear();
                this->outer.clear();
                this->values.clear();
            }

            /**
             * @brief Returns the compressed state.
             *
             * @return true
             * @return false
             */
            inline bool is_compressed() const {
                return this->compressed;
            }

            // OPERATIONS.

            /**
             * @brief Sparse matrix unary +.
             * 
             * @return Sparse 
             */
            Sparse operator +() const {
                return *this;
            }

            /**
             * @brief Sparse matrix unary -.
             * 
             * @return Sparse 
             */
            Sparse operator -() const {
                Sparse result{*this};

                if(!(result.compressed))
                    for(auto &[key, element]: result.elements)
                        element = -element;
                else
                    for(auto &value: result.values)
                        value = -value;
                
                return result;
            }

            /**
             * @brief Sparse matrix sum.
             * 
             * @param sparse 
             * @return Sparse 
             */
            Sparse operator +(const Sparse &sparse) const {
                #ifndef NDEBUG // Integrity checks.
                assert((this->rows == sparse.rows) && (this->columns == sparse.columns));
                #endif

                Sparse result{*this};
                result.uncompress();

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        result.add(key[0], key[1], element);
                else
                    for(std::size_t j = 0; j < sparse.rows; ++j)
                        for(std::size_t k = sparse.inner[j]; k < sparse.inner[j + 1]; ++k)
                            result.add(j, sparse.outer[k], sparse.values[k]);

                return result;
            }

            /**
             * @brief Sparse matrix sum and assignation.
             * 
             * @param sparse 
             * @return Sparse& 
             */
            Sparse &operator +=(const Sparse &sparse) {
                #ifndef NDEBUG // Integrity checks.
                assert((this->rows == sparse.rows) && (this->columns == sparse.columns));
                assert(!this->compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        this->add(key[0], key[1], element);
                else
                    for(std::size_t j = 0; j < sparse.rows; ++j)
                        for(std::size_t k = sparse.inner[j]; k < sparse.inner[j + 1]; ++k)
                            this->add(j, sparse.outer[k], sparse.values[k]);

                return *this;
            }

            /**
             * @brief Sparse matrix difference.
             * 
             * @param sparse 
             * @return Sparse 
             */
            Sparse operator -(const Sparse &sparse) const {
                #ifndef NDEBUG // Integrity checks.
                assert((this->rows == sparse.rows) && (this->columns == sparse.columns));
                #endif

                Sparse result{*this};
                result.uncompress();

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        result.add(key[0], key[1], -element);
                else
                    for(std::size_t j = 0; j < sparse.rows; ++j)
                        for(std::size_t k = sparse.inner[j]; k < sparse.inner[j + 1]; ++k)
                            result.add(j, sparse.outer[k], -sparse.values[k]);

                return result;
            }

            /**
             * @brief Sparse matrix difference and assignation.
             * 
             * @param sparse 
             * @return Sparse& 
             */
            Sparse &operator -=(const Sparse &sparse) {
                #ifndef NDEBUG // Integrity checks.
                assert((this->rows == sparse.rows) && (this->columns == sparse.columns));
                assert(!this->compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        this->add(key[0], key[1], -element);
                else
                    for(std::size_t j = 0; j < sparse.rows; ++j)
                        for(std::size_t k = sparse.inner[j]; k < sparse.inner[j + 1]; ++k)
                            this->add(j, sparse.outer[k], -sparse.values[k]);

                return *this;
            }

            /**
             * @brief Sparse matrix scalar product.
             * 
             * @param scalar 
             * @return Sparse 
             */
            Sparse operator *(const T &scalar) const {
                Sparse result{*this};

                if(!(this->compressed)) 
                    for(auto &[key, element]: result.elements)
                        element *= scalar;
                else {
                    
                    for(auto &value: result.values)
                        value *= scalar;

                }

                return result;
            }

            /**
             * @brief Friend Sparse matrix scalar product.
             * 
             * @param scalar 
             * @param sparse 
             * @return Sparse 
             */
            friend Sparse operator *(const T &scalar, const Sparse &sparse) {
                Sparse result{sparse};

                if(!(sparse.compressed))
                    for(auto &[key, element]: result.elements)
                        element *= scalar;
                else {
                    
                    for(auto &value: result.values)
                        value *= scalar;
                }

                return result;
            }

            /**
             * @brief Sparse matrix scalar product and assignation.
             * 
             * @param scalar 
             * @return Sparse& 
             */
            Sparse &operator *=(const T &scalar) {
                if(!(this->compressed))
                    for(auto &[key, element]: this->elements)
                        element *= scalar;
                else {
                    
                    for(std::size_t j = 0; j < this->values.size(); ++j)
                        this->values[j] *= scalar;
                }

                return *this;
            }

            /**
             * @brief Sparse matrix * Vector product.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> operator *(const Vector<T> &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->columns == vector.length);
                #endif

                Vector<T> result{this->rows};

                if(!(this->compressed))
                    for(const auto &[key, element]: this->elements)
                        result[key[0]] += element * vector[key[1]];
                else
                    for(std::size_t j = 0; j < this->rows; ++j)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            result[j] += this->values[k] * vector[this->outer[k]];

                return result;
            }

            /**
             * @brief Sparse matrix * Sparse matrix.
             * 
             * @param sparse 
             * @return Sparse 
             */
            Sparse operator *(const Sparse &sparse) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->columns == sparse.rows);
                #endif

                Sparse result{this->rows, sparse.columns};

                // Temporary solution, may be extremely slow.
                for(std::size_t j = 0; j < this->rows; ++j)
                    for(std::size_t k = 0; k < sparse.columns; ++k)
                        result.insert(j, k, dot(this->row(j), sparse.column(k)));

                return result;
            }

            // LINEAR.

            /**
             * @brief Solves a linear system in the form of Ax = b.
             * 
             * @tparam S 
             * @param vector 
             * @return Vector<T> 
             */
            template<Solver S = Conjugate>
            Vector<T> solve(const Vector<T> &vector) const {
                
                if constexpr (S == Conjugate)
                    return this->conjugate_gradient(vector);

                if constexpr (S == Descent)
                    return this->gradient_descent(vector);

                if constexpr (S == Gauss)
                    return this->gauss_seidel(vector);

                if constexpr (S == Kaczmarz)
                    return this->kaczmarz(vector);

                if constexpr (S == RandomKaczmarz)
                    return this->randomized_kaczmarz(vector);

                // Default.
                return this->conjugate_gradient(vector);
            }

            // OUTPUT.

            /**
             * @brief Sparse matrix output.
             * 
             * @param ost 
             * @param sparse 
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Sparse &sparse) {
                if(!(sparse.compressed))
                    for(const auto &[key, element]: sparse.elements) {
                        ost << "(" << key[0] << ", " << key[1] << "): " << element;

                        if(key != (*--sparse.elements.end()).first)
                            ost << std::endl;
                    }
                else {
                    for(std::size_t j = 0; j < sparse.rows; ++j)
                        for(std::size_t k = sparse.inner[j]; k < sparse.inner[j + 1]; ++k) {
                            ost << "(" << j << ", " << sparse.outer[k] << "): " << sparse.values[k];

                            if(k < sparse.inner[sparse.rows] - 1)
                                ost << std::endl;
                        }
                }

                return ost;
            }

        private:

            // SOLVERS.

            /**
             * @brief Conjugate Gradient method.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> conjugate_gradient(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                assert(this->columns == vector.length);
                assert(std::abs(mtrace(*this)) > TOLERANCE);
                #endif

                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Problem's size.
                const std::size_t size = this->rows;
                
                // Iterations.
                std::size_t iterations = 0;

                // Solution.
                Vector<T> solution{size}, old_solution{size};

                // Search.
                Vector<T> search{size};

                // Parameters.
                T alpha, beta;

                // Direction.
                Vector<T> direction = vector;

                // Residual.
                Vector<T> residual = vector, old_residual = vector;

                // Target.
                Sparse target{*this};
                target.compress();

                // Method.
                do {
                    ++iterations;

                    #ifndef NVERBOSE
                    if(!(iterations % 50))
                        std::cout << "\tConjugate Gradient, iteration: " << iterations << std::endl;
                    #endif

                    // Search.
                    search = target * direction;

                    // Old solution and residual.
                    old_solution = solution;
                    old_residual = residual;

                    // Alpha computation.
                    alpha = dot(residual, residual) / dot(direction, search);

                    // Step.
                    solution += alpha * direction;
                    residual -= alpha * search;

                    // Checks.
                    if(((vector - target * solution).norm() <= ALGEBRA_TOLERANCE) || ((old_solution - solution).norm() <= ALGEBRA_TOLERANCE))
                        break;

                    //  Beta computation.
                    beta = dot(residual, residual) / dot(old_residual, old_residual);

                    // Direction evaluation.
                    direction = residual + beta * direction;

                } while(iterations < ALGEBRA_ITER_MAX);

                #ifndef NVERBOSE
                std::cout << "Results:" << std::endl;
                std::cout << "\tIterations: " << iterations << std::endl;
                std::cout << "\tResidual: " << residual.norm() << std::endl;
                #endif

                return solution;
            }

            /**
             * @brief Gradient Descent method.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> gradient_descent(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                assert(this->columns == vector.length);
                assert(std::abs(mtrace(*this)) > TOLERANCE);
                #endif

                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Problem's size.
                const std::size_t size = this->rows;
                
                // Iterations.
                std::size_t iterations = 0;

                // Solution.
                Vector<T> solution{size}, old_solution{size};

                // Parameters.
                T alpha;

                // Residual.
                Vector<T> residual = vector, old_residual = vector;

                // Target.
                Sparse target{*this};
                target.compress();

                // Method.
                do {
                    ++iterations;

                    #ifndef NVERBOSE
                    if(!(iterations % 50))
                        std::cout << "\tGradient Descent, iteration: " << iterations << std::endl;
                    #endif

                    old_solution = solution;

                    // Alpha computation.
                    alpha = dot(residual, residual) / dot(residual, target * residual);

                    // Step.
                    solution += alpha * residual;
                    residual = vector - target * solution;

                    // Checks.
                    if((residual.norm() <= ALGEBRA_TOLERANCE) || ((old_solution - solution).norm() <= ALGEBRA_TOLERANCE))
                        break;

                    residual = vector - target * solution;

                } while(iterations < ALGEBRA_ITER_MAX);

                #ifndef NVERBOSE
                std::cout << "Results:" << std::endl;
                std::cout << "\tIterations: " << iterations << std::endl;
                std::cout << "\tResidual: " << residual.norm() << std::endl;
                #endif

                return solution;
            }

            /**
             * @brief Gauss-Seidel method.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> gauss_seidel(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                assert(this->columns == vector.length);
                assert(std::abs(mtrace(*this)) > TOLERANCE);
                #endif

                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Problem's size.
                const std::size_t size = this->rows;

                // Target.
                Sparse target{*this};
                target.compress();
                
                // Iterations.
                std::size_t iterations = 0;

                // Splitting.
                Sparse lower = this->lower() + this->diagonal();
                Sparse upper = this->upper();

                // Inversion of lower.
                Sparse lower_inv{size, size};
                lower.compress();

                for(std::size_t h = 0; h < size; ++h) {
                    Vector<T> column{size};

                    for(std::size_t j = 0; j < size; ++j) {
                        for(std::size_t k = lower.inner[j]; k < lower.inner[j + 1]; ++k) {
                            T sum = static_cast<T>(0);

                            for(std::size_t l = lower.inner[j]; l < lower.inner[j + 1] - 1; ++l)
                                sum += lower.values[l] * column[lower.outer[l]];

                            column[j] = (static_cast<T>(j == h) - sum) / lower.values[lower.inner[j + 1] - 1];
                        }
                    }

                    for(std::size_t k = 0; k < size; ++k) // May introduce fill-ins.
                        lower_inv.elements[{k, h}] = column[k];
                }

                // T-matrix and C-vector.
                Sparse t_matrix = - (lower_inv * upper);
                Vector<T> c_vector = lower_inv * vector;
                t_matrix.compress();

                // Solution.
                Vector<T> solution{size}, old_solution{size};

                // Method.
                do {
                    ++iterations;

                    #ifndef NVERBOSE
                    if(!(iterations % 50))
                        std::cout << "\tGauss-Seidel, iteration: " << iterations << std::endl;
                    #endif

                    // Step.
                    old_solution = solution;
                    solution = t_matrix * solution + c_vector;

                    // Checks.
                    if((target * solution - vector).norm() <= ALGEBRA_TOLERANCE)
                        break;

                } while(iterations < ALGEBRA_ITER_MAX);

                #ifndef NVERBOSE
                std::cout << "Results:" << std::endl;
                std::cout << "\tIterations: " << iterations << std::endl;
                std::cout << "\tResidual: " << (target * solution - vector).norm() << std::endl;
                #endif

                return solution;
                
            }

            /**
             * @brief Kaczmarz method.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> kaczmarz(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                assert(this->columns == vector.length);
                assert(std::abs(mtrace(*this)) > TOLERANCE);
                #endif

                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Problem's size.
                const std::size_t size = this->rows;
                
                // Iterations.
                std::size_t iterations = 0;

                // Solution.
                Vector<T> solution{size}, old_solution{size};

                // Target.
                Sparse target{*this};
                target.compress();

                // Method.
                do {
                    ++iterations;

                    #ifndef NVERBOSE
                    if(!(iterations % 50))
                        std::cout << "\tKaczmarz, iteration: " << iterations << std::endl;
                    #endif

                    // Index.
                    std::size_t index = iterations % size;

                    // Step.
                    Vector<T> row{size};
                    T product = static_cast<T>(0);
                    Real norm = 0.0;

                    for(std::size_t k = target.inner[index]; k < target.inner[index + 1]; ++k) {
                        const T value = target.values[k];

                        row[target.outer[k]] = value;
                        norm += std::abs(value) * std::abs(value);
                        product += solution[target.outer[k]] * value;
                    }

                    old_solution = solution;
                    solution += ((vector[index] - product) / norm) * row;

                    // Checks.
                    if((target * solution - vector).norm() <= ALGEBRA_TOLERANCE)
                        break;

                } while(iterations < ALGEBRA_ITER_MAX);

                #ifndef NVERBOSE
                std::cout << "Results:" << std::endl;
                std::cout << "\tIterations: " << iterations << std::endl;
                std::cout << "\tResidual: " << (target * solution - vector).norm() << std::endl;
                #endif

                return solution;
            }

            /**
             * @brief Randomized Kaczmarz method.
             * 
             * @param vector 
             * @return Vector<T> 
             */
            Vector<T> randomized_kaczmarz(const Vector<T> &vector) const {
                #ifndef NDEBUG
                assert(this->rows == this->columns);
                assert(this->columns == vector.length);
                assert(std::abs(mtrace(*this)) > TOLERANCE);
                #endif

                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Problem's size.
                const std::size_t size = this->rows;
                
                // Iterations.
                std::size_t iterations = 0;

                // Solution.
                Vector<T> solution{size}, old_solution{size};

                // Target.
                Sparse target{*this};
                target.compress();

                // Probabilities.
                Vector<Real> norms{size}, probabilities{size + 1};
                Real norm = 0.0;

                for(std::size_t j = 0; j < size; ++j) {
                    for(std::size_t k = target.inner[j]; k < target.inner[j + 1]; ++k)
                        norms[j] += std::abs(target.values[k]) * std::abs(target.values[k]);

                    norm += norms[j];
                }

                for(std::size_t j = 0; j < size + 1; ++j) {
                    Real sum = 0.0;

                    for(std::size_t k = 0; k < j; ++k)
                        sum += norms[k] / norm;

                    probabilities[j] = sum;
                }
                        
                #ifndef NVERBOSE
                std::cout << "Solving a linear system." << std::endl;
                #endif

                // Method.
                do {
                    ++iterations;

                    #ifndef NVERBOSE
                    if(!(iterations % 50))
                        std::cout << "\tRandomized Kaczmarz, iteration: " << iterations << std::endl;
                    #endif

                    // Index.
                    std::size_t index = iterations % size; // Fallback.
                    Real random = static_cast<Real>(std::rand()) / static_cast<Real>(RAND_MAX);

                    for(std::size_t j = 0; j < size; ++j)
                        if((probabilities[j] <= random) && (random <= probabilities[j + 1])) {
                            index = j;
                            break;
                        }

                    // Step.
                    Vector<T> row{size};
                    T product = static_cast<T>(0);
                    Real norm = 0.0;

                    for(std::size_t k = target.inner[index]; k < target.inner[index + 1]; ++k) {
                        const T value = target.values[k];

                        row[target.outer[k]] = value;
                        norm += std::abs(value) * std::abs(value);
                        product += solution[target.outer[k]] * value;
                    }

                    old_solution = solution;
                    solution += ((vector[index] - product) / norm) * row;

                    // Checks.
                    if((target * solution - vector).norm() <= ALGEBRA_TOLERANCE)
                        break;

                } while(iterations < ALGEBRA_ITER_MAX);

                #ifndef NVERBOSE
                std::cout << "Results:" << std::endl;
                std::cout << "\tIterations: " << iterations << std::endl;
                std::cout << "\tResidual: " << (target * solution - vector).norm() << std::endl;
                #endif

                return solution;
            }
    };

    /**
     * @brief Multiplicative trace.
     * 
     * @tparam T 
     * @param matrix 
     * @return T 
     */
    template<NumericType T>
    T mtrace(const Sparse<T> &sparse) {
        #ifndef NDEBUG // Integrity check.
        assert(sparse.rows == sparse.columns);
        #endif

        T product = static_cast<T>(1);

        for(std::size_t j = 0; j < sparse.rows; ++j)
            product *= sparse(j, j); // Slow.

        return product;
    }

}

#endif