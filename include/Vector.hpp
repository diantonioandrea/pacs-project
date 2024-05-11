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

// Copy.
#include <algorithm>

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
            Vector(const std::size_t &length, const std::vector<T> &elements): length{length} {
                #ifndef NDEBUG // Integrity check.
                assert(length > 0);
                assert(elements.size() == length);
                #endif

                this->elements.resize(elements.size());
                std::ranges::copy(elements.begin(), elements.end(), this->elements.begin());
            }

            /**
             * @brief Copy constructor.
             * 
             * @param vector 
             */
            Vector(const Vector &vector): length{vector.length} {
                this->elements.resize(vector.length);
                std::ranges::copy(vector.elements.begin(), vector.elements.end(), this->elements.begin());
            }
            
            /**
             * @brief Copy operator.
             * 
             * @param vector 
             * @return Vector& 
             */
            Vector &operator =(const Vector &vector) {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                this->elements.resize(vector.length);
                std::ranges::copy(vector.elements.begin(), vector.elements.end(), this->elements.begin());

                return *this;
            }

            /**
             * @brief Scalar copy operator.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator =(const T &scalar) {
                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] = scalar;

                return *this;
            }

            // CONVERSION.

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

            /**
             * @brief Returns the [j, k) range.
             * 
             * @param j 
             * @param k 
             * @return Vector 
             */
            Vector operator ()(const std::size_t &j, const std::size_t &k) const {
                #ifndef NDEBUG // Integrity check.
                assert(j != k);
                assert((j < this->length) && (k < this->length + 1));
                #endif

                Vector result{(j < k) ? k - j : j - k};

                if(j < k) {
                    for(std::size_t h = j; h < k; ++h)
                        result[h - j] = this->elements[h];
                } else {
                    for(std::size_t h = k; h > j; h--)
                        result[h - k] = this->elements[h];
                }

                return result;
            }

            /**
             * @brief Returns the [j, end) range.
             * 
             * @param j 
             * @return Vector 
             */
            Vector operator ()(const std::size_t &j) const {
                #ifndef NDEBUG // Integrity check.
                assert(j < this->length);
                #endif

                Vector result{this->length - j};
                
                for(std::size_t h = j; h < this->length; ++h)
                    result[h - j] = this->elements[h];

                return result;
            }
            
            /**
             * @brief Returns a sub-Vector given a std::vector of indices.
             * 
             * @param indices 
             * @return Vector 
             */
            Vector operator ()(const std::vector<std::size_t> &indices) const {
                #ifndef NDEBUG // Integrity check.
                for(const auto &index: indices)
                    assert(index < this->length);
                #endif

                Vector result{indices.size()};

                for(std::size_t j = 0; j < indices.size(); ++j)
                    result[j] = this->elements[indices[j]];

                return result;
            }

            /**
             * @brief Sets a sub-Vector given a std::vector of indices and a Vector of values.
             * 
             * @param indices 
             * @param values 
             * @return Vector 
             */
            void operator ()(const std::vector<std::size_t> &indices, const Vector &values) {
                #ifndef NDEBUG // Integrity check.
                assert(indices.size() == values.length);
                for(const auto &index: indices)
                    assert(index < this->length);
                #endif

                for(std::size_t j = 0; j < indices.size(); ++j)
                    this->elements[indices[j]] = values[j];
            }

            // OPERATIONS.

            /**
             * @brief Vector unary +.
             * 
             * @return Vector 
             */
            Vector operator +() const {
                return *this;
            }

            /**
             * @brief Vector unary -.
             * 
             * @return Vector 
             */
            Vector operator -() const {
                Vector result{*this};

                for(auto &element: result.elements)
                    element = -element;

                return result;
            }

            /**
             * @brief Scalar product.
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator *(const T &scalar) const {
                Vector result{*this};

                for(auto &element: result.elements)
                    element *= scalar;

                return result;
            }

            /**
             * @brief Friend scalar product.
             * 
             * @param scalar 
             * @param vector 
             * @return Vector 
             */
            friend Vector operator *(const T &scalar, const Vector &vector) {
                Vector result{vector};

                for(auto &element: result.elements)
                    element *= scalar;

                return result;
            }

            /**
             * @brief Scalar product and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator *=(const T &scalar) {
                for(auto &element: this->elements)
                    element *= scalar;

                return *this;
            }

            /**
             * @brief Scalar division.
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator /(const T &scalar) const {
                Vector result{*this};

                for(auto &element: result.elements)
                    element /= scalar;

                return result;
            }

            /**
             * @brief Friend scalar division.
             * 
             * @param scalar 
             * @param vector 
             * @return Vector 
             */
            friend Vector operator /(const T &scalar, const Vector &vector) {
                Vector result{vector.length, scalar};

                for(std::size_t j = 0; j < result.length; ++j)
                    result.elements[j] /= vector.elements[j];

                return result;
            }

            /**
             * @brief Scalar division and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator /=(const T &scalar) {
                for(auto &element: this->elements)
                    element /= scalar;

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

                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] += vector.elements[j];

                return *this;
            }

            /**
             * @brief Scalar "sum".
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator +(const T &scalar) const {
                Vector result{*this};

                for(auto &element: result.elements)
                    element += scalar;
                
                return result;
            }

            /**
             * @brief Friend scalar "sum".
             * 
             * @param scalar 
             * @param vector 
             * @return Vector 
             */
            friend Vector operator +(const T &scalar, const Vector &vector) {
                Vector result{vector};

                for(auto &element: result.elements)
                    element += scalar;
                
                return result;
            }

            /**
             * @brief Scalar "sum" and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator +=(const T &scalar) {
                for(auto &element: this->elements)
                    element += scalar;
                
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

                for(std::size_t j = 0; j < this->length; ++j)
                    this->elements[j] -= vector.elements[j];

                return *this;
            }

            /**
             * @brief Scalar "difference".
             * 
             * @param scalar 
             * @return Vector 
             */
            Vector operator -(const T &scalar) const {
                Vector result{*this};

                for(auto &element: result.elements)
                    element -= scalar;
                
                return result;
            }

            /**
             * @brief Friend scalar "difference".
             * 
             * @param scalar 
             * @param vector 
             * @return Vector 
             */
            friend Vector operator -(const T &scalar, const Vector &vector) {
                Vector result{vector};

                for(auto &element: result.elements)
                    element = scalar - element;
                
                return result;
            }

            /**
             * @brief Scalar "difference" and assignation.
             * 
             * @param scalar 
             * @return Vector& 
             */
            Vector &operator -=(const T &scalar) {
                for(auto &element: this->elements)
                    element -= scalar;
                
                return *this;
            }

            /**
             * @brief Vector element-wise product.
             * 
             * @param vector 
             * @return T 
             */
            Vector operator *(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                Vector result{*this};

                for(std::size_t j = 0; j < result.length; ++j)
                    result.elements[j] *= vector.elements[j];

                return result;
            }

            /**
             * @brief Vector element-wise division.
             * 
             * @param vector 
             * @return T 
             */
            Vector operator /(const Vector &vector) const {
                #ifndef NDEBUG // Integrity check.
                assert(this->length == vector.length);
                #endif

                Vector result{*this};

                for(std::size_t j = 0; j < result.length; ++j)
                    result.elements[j] /= vector.elements[j];

                return result;
            }

            // NORM.

            /**
             * @brief Returns the l2 norm of the Vector.
             * 
             * @return Real 
             */
            Real norm() const {
                Real norm = 0.0;

                for(const auto &element: this->elements)
                    norm += std::abs(element) * std::abs(element);

                return std::sqrt(norm);
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

    /**
     * @brief Dot product.
     * 
     * @tparam T 
     * @param first 
     * @param second 
     * @return T 
     */
    template<NumericType T>
    T dot(const Vector<T> &first, const Vector<T> &second) {
        #ifndef NDEBUG
        assert(first.length == second.length);
        #endif
    
        T result = static_cast<T>(0);

        for(std::size_t j = 0; j < first.length; ++j)
            result += first[j] * second[j];

        return result;
    }

    /**
     * @brief Stacks two vectors.
     * 
     * @tparam T 
     * @param first 
     * @param second 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> stack(const Vector<T> &first, const Vector<T> &second) {
        Vector<T> result{first.length + second.length};

        {
            for(std::size_t j = 0; j < first.length; ++j)
                result[j] = first[j];

            for(std::size_t j = first.length; j < first.length + second.length; ++j)
                result[j] = second[j - first.length];
        }

        return result;
    }

    /**
     * @brief Flips a vector.
     * 
     * @tparam T 
     * @param vector 
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> flip(const Vector<T> &vector) {
        Vector<T> result{vector.length};

        for(std::size_t j = 0; j < result.length; ++j)
            result[j] = vector[vector.length - j - 1];

        return result;
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

        for(std::size_t j = 0; j < vector.length; ++j)
            result[j] = std::abs(vector[j]);

        return result;
    }

}

#endif