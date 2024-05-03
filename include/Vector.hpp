/**
 * @file Vector.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef VECTOR_PACS
#define VECTOR_PACS

// Type.
#include <Type.hpp>

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

namespace pacs {

    template<NumericType T>
    class Vector {
        protected:

            // Elements.
            std::vector<T> elements;

        public:

            // Shape.
            const std::size_t length;

            // CONSTRUCTORS.

            /**
             * @brief Constructs a new empty Vector.
             * 
             * @param length 
             */
            Vector(const std::size_t &length): length{length} {
                #ifndef NDEBUG // Integrity check.
                assert(length > 0);
                #endif

                this->elements.resize(length, static_cast<T>(0));
            }

            /**
             * @brief Constructs a new Vector from a given std::vector.
             * 
             * @param length 
             * @param elements 
             */
            Vector(const std::size_t &length, const std::vector<T> &elements): length{length}, elements{elements} {
                #ifndef NDEBUG // Integrity check.
                assert(length > 0);
                assert(elements.size() == length);
                #endif
            }

            /**
             * @brief Copy constructor.
             * 
             * @param vector 
             */
            Vector(const Vector &vector): length{vector.length}, elements{vector.elements} {}
            
            Vector &operator =(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                this->elements = vector.elements;
            }

            // CONVERSION

            /**
             * @brief Converts the Vector into a std::vector<T>.
             * 
             * @return std::vector<T> 
             */
            operator std::vector<T>() const {
                return this->elements;
            }

            // OPERATORS.

            /**
             * @brief Vector scalar product.
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator *(const T &scalar) const {
                Vector result{*this};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    result.elements[j] *= scalar;

                return result;
            }

            /**
             * @brief Vector scalar product and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator *(const T &scalar) {
                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] *= scalar;

                return *this;
            }

            /**
             * @brief Vector sum.
             * 
             * @param vector 
             * @return Vector 
             */
            Vector operator +(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                Vector result{*this};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    result.elements[j] += vector.elements[j];

                return result;
            }

            /**
             * @brief Vector sum and assignation.
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator +=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] += vector.elements[j];

                return *this;
            }

            /**
             * @brief Vector difference.
             * 
             * @param vector 
             * @return Vector 
             */
            Vector operator -(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                Vector result{*this};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    result.elements[j] -= vector.elements[j];

                return result;
            }

            /**
             * @brief Vector difference and assignation.
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator -=(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] -= vector.elements[j];

                return *this;
            }

            T operator *(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                T result = static_cast<T>(0);

                #pragma omp parallel for reduction(+: result)
                for(std::size_t j = 0; j < this->length; ++j)
                    result += this->elements[j] * vector.elements[j];

                return result;
            }

            // OUTPUT.

            /**
             * @brief Vector output.
             * 
             * @param ost 
             * @param vector 
             * @return std::ostream& 
             */
            friend std::ostream &operator <<(std::ostream &ost, const Vector &vector) {
                for(const auto &element: vector.elements)
                    ost << element << " ";

                return ost;
            }
    };

}

#endif