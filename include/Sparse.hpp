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

namespace pacs {
    
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

                #pragma omp parallel for
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

                #pragma omp parallel for
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
                        for(std::size_t h = this->inner[j]; h < this->inner[j + 1]; ++j)
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

                #pragma omp parallel for
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

                #pragma omp parallel for
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

                        #ifndef NDEBUG
                        if(std::abs(value) > TOLERANCE) {
                            this->outer.emplace_back(key[1]);
                            this->values.emplace_back(value);
                            ++index;
                        }
                        #else
                        this->outer.emplace_back(key[1]);
                        this->values.emplace_back(value);
                        ++index;
                        #endif

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
                
                #ifdef DYNAMIC_SPARSE
                result.uncompress();
                #endif

                #ifndef NDEBUG // Integrity check.
                assert(!result.compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        result.add(key[0], key[1], element);
                else
                    #pragma omp parallel for collapse(2)
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
                #endif

                #ifdef DYNAMIC_SPARSE
                this->uncompress();
                #endif

                #ifndef NDEBUG // Integrity check.
                assert(!this->compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        this->add(key[0], key[1], element);
                else
                    #pragma omp parallel for collapse(2)
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
                
                #ifdef DYNAMIC_SPARSE
                result.uncompress();
                #endif

                #ifndef NDEBUG // Integrity check.
                assert(!result.compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        result.add(key[0], key[1], -element);
                else
                    #pragma omp parallel for collapse(2)
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
                #endif

                #ifdef DYNAMIC_SPARSE
                this->uncompress();
                #endif

                #ifndef NDEBUG // Integrity check.
                assert(!this->compressed);
                #endif

                if(!(sparse.compressed))
                    for(auto &[key, element]: sparse.elements)
                        this->add(key[0], key[1], -element);
                else
                    #pragma omp parallel for collapse(2)
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
                    
                    #pragma omp parallel for
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
                    
                    #pragma omp parallel for
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
                    
                    #pragma omp parallel for
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
                else {
                    for(std::size_t j = 0; j < this->rows; ++j) {
                        T product = static_cast<T>(0);

                        #pragma omp parallel for reduction(+: product)
                        for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                            product += this->values[k] * vector[this->outer[k]];

                        result[j] = product;
                    }
                }

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
                #pragma omp parallel for collapse(2)
                for(std::size_t j = 0; j < this->rows; ++j)
                    for(std::size_t k = 0; k < sparse.columns; ++k)
                        result.insert(j, k, dot(this->row(j), sparse.column(k)));

                #ifdef DYNAMIC_SPARSE
                if(this->compressed && sparse.compressed)
                    result.compress();
                #endif

                return result;
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

        #pragma omp parallel for reduction(*: product)
        for(std::size_t j = 0; j < sparse.rows; ++j)
            product *= sparse(j, j); // Slow.

        return product;
    }

}

#endif