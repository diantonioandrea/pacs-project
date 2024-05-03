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
            Sparse(const Sparse &sparse): rows{sparse.rows}, columns{sparse.columns}, compressed{sparse.compressed} {
                if(!(sparse.compressed)) {
                    this->elements = sparse.elements;
                } else {
                    this->inner = sparse.inner;
                    this->outer = sparse.outer;
                    this->values = sparse.values;
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

                if(!(sparse.compressed)) {
                    this->elements = sparse.elements;

                    this->inner.clear();
                    this->outer.clear();
                    this->values.clear();
                } else {
                    this->elements.clear();

                    this->inner = sparse.inner;
                    this->outer = sparse.outer;
                    this->values = sparse.values;
                }

                return *this;
            }

            // READ AND WRITE.

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
                for(std::size_t i = this->inner[j]; i < this->inner[j + 1]; ++i) {
                    if(k == this->outer[i])
                        return this->values[i];
                }

                // Default return.
                return static_cast<T>(0);
            }

            /**
             * @brief Insert a new element.
             * 
             * @param j 
             * @param k 
             * @param element 
             */
            void insert(const std::size_t &j, const std::size_t &k, const T &element) {
                #ifndef NDEBUG // Out-of-bound and uncompression check.
                assert((j < rpws) && (k < columns));
                assert(!(this->compressed));
                
                if(std::abs(element) > TOLERANCE)
                    this->elements[{j, k}] = element;
                #else
                this->elements[{j, k}] = element;
                #endif
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
                    for(auto it = this->elements.lower_bound(current); (*it).rows < (*(this->elements.lower_bound(next))).rows; ++it) {
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
                for(std::size_t j = 0; j < this->inner.size() - 1; ++j) {
                    for(std::size_t k = this->inner[j]; k < this->inner[j + 1]; ++k)
                        this->elements[{j, this->outer[k]}] = this->values[k];
                }

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

            // OPERATORS

            // ...

            // OUTPUT.

            /**
             * @brief Sparse matrix output.
             * 
             * @param ost 
             * @param sparse 
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Sparse &sparse) {
                if(!(sparse.compressed)) {
                    for(const auto &[key, value]: sparse.elements) {
                        ost << "(" << key[0] << ", " << key[1] << "): " << value;

                        if(key != (*--sparse.elements.end()).first)
                            ost << std::endl;
                    }

                } else {
                    ost << "Inner: ";

                    for(const auto &value: sparse.inner)
                        ost << value << " ";

                    ost << std::endl << "Outer: ";

                    for(const auto &value: sparse.outer)
                        ost << value << " ";

                    ost << std::endl <<  "Values: ";

                    for(const auto &value: sparse.values)
                        ost << value << " ";
                }

                return ost;
            }
    };

}

#endif