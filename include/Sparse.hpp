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

namespace pacs {
    
    /**
     * @brief Sparse matrix class.
     * 
     * @tparam T Matrix' type.
     */
    template<NumericType T>
    class Sparse {
        private:

            // Shape.
            const std::size_t rows;
            const std::size_t columns;

            // Compression flag.
            bool compressed = false;

            // COOmap dynamic storage format.
            mutable std::map<std::array<std::size_t, 2>, T> elements;

            // CSR compressed storage format.
            std::vector<std::size_t> inner;
            std::vector<std::size_t> outer;
            std::vector<T> values;

        public:

            // CONSTRUCTORS.

            /**
             * @brief Construct a new empty Sparse matrix.
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
             * @brief Construct a new Sparse matrix from a given std::map.
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
             * @brief Construct a new Sparse matrix from given inner, outer and values vectors.
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
    };

}

#endif