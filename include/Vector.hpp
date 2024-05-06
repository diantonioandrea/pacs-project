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
             * @brief Constructs a new homogeneous Vector.
             * 
             * @param length 
             * @param value 
             */
            Vector(const std::size_t &length, const T &value): length{length} {
                #ifndef NDEBUG // Integrity check.
                assert(length > 0);
                #endif

                this->elements.resize(length, value);
            }

            /**
             * @brief Constructs a new Vector from a given std::vector.
             * 
             * @param length 
             * @param elements 
             */
            Vector(const std::size_t &length, const std::vector<T> &elements): elements{elements}, length{length} {
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
            Vector(const Vector &vector): elements{vector.elements}, length{vector.length} {}
            
            Vector &operator =(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                this->elements = vector.elements;

                return *this;
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

            // READ AND WRITE.

            /**
             * @brief Const subscript operator, returns the j-th element.
             * 
             * @param j 
             * @return T 
             */
            inline T operator [](const std::size_t &j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->length);
                #endif

                return this->elements[j];
            }

            /**
             * @brief Subscript operator, returns a reference to the j-th element.
             * 
             * @param j 
             * @return T& 
             */
            inline T &operator [](const std::size_t &j) {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->length);
                #endif

                return this->elements[j];
            }

            // OPERATORS.

            /**
             * @brief Vector scalar product.
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator *(const T &scalar) const {
                Vector result{this->length};

                #pragma omp parallel for
                for(std::size_t j = 0; j < this->length; ++j)
                    result[j] = this->elements[j] * scalar;

                return result;
            }

            /**
             * @brief Friend Vector scalar product.
             * 
             * @param scalar 
             * @param vector 
             * @return Vector 
             */
            friend Vector operator *(const T &scalar, const Vector &vector) {
                Vector result{vector.length};

                #pragma omp parallel for
                for(std::size_t j = 0; j < vector.length; ++j)
                    result[j] = vector.elements[j] * scalar;

                return result;
            }

            /**
             * @brief Vector scalar product and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator *=(const T &scalar) {
                #pragma omp parallel for
                for(auto &element: this->elements)
                    element *= scalar;

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

            /**
             * @brief Vector dot product.
             * 
             * @param vector 
             * @return T 
             */
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

    // METHODS.

    /**
     * @brief Maximum value of a Vector.
     * 
     * @tparam T 
     * @param vector 
     * @return T 
     */
    template<NumericType T>
    T max(const Vector<T> &vector) {
        T max = vector[0];

        for(std::size_t j = 1; j < vector.length; ++j)
            max = (vector[j] > max) ? vector[j] : max;

        return max;
    }

    /**
     * @brief Minimum value of a Vector.
     * 
     * @tparam T 
     * @param vector 
     * @return T 
     */
    template<NumericType T>
    T min(const Vector<T> &vector) {
        T min = vector[0];

        for(std::size_t j = 1; j < vector.length; ++j)
            min = (vector[j] < min) ? vector[j] : min;

        return min;
    }

}

namespace std {

    /**
     * @brief std::abs overload for Vectors.
     * 
     * @tparam T 
     * @param vector 
     * @return pacs::Vector<T> 
     */
    template<pacs::NumericType T>
    pacs::Vector<T> abs(const pacs::Vector<T> vector) {
        pacs::Vector<T> result{vector.length};

        #pragma omp parallel for
        for(std::size_t j = 0; j < vector.length; ++j)
            result[j] = std::abs(vector[j]);

        return result;
    }

}

#endif