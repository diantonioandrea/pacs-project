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

#include <Type.hpp>

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace pacs {

    /**
     * @brief Vector structure.
     * 
     * @tparam T 
     */
    template<NumericType T>
    struct Vector {

        // Shape.
        const std::size_t length;

        // Elements.
        std::vector<T> elements;

        // CONSTRUCTORS.

        /**
         * @brief Constructs a new empty Vector.
         * 
         * @param length 
         */
        Vector(const std::size_t &length): length{length}, elements(length, static_cast<T>(0)) {
            #ifndef NDEBUG // Integrity check.
            assert(length > 0);
            #endif
        }

        /**
         * @brief Constructs a new homogeneous Vector.
         * 
         * @param length 
         * @param value 
         */
        Vector(const std::size_t &length, const T &value): length{length}, elements(length, value) {
            #ifndef NDEBUG // Integrity check.
            assert(length > 0);
            #endif
        }

        /**
         * @brief Constructs a new Vector from a given std::vector.
         * 
         * @param length 
         * @param elements 
         */
        Vector(const std::size_t &length, const std::vector<T> &elements): length{length}, elements(elements.begin(), elements.end()) {
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
        Vector(const Vector &vector): length{vector.length}, elements(vector.elements.begin(), vector.elements.end()) {}
        
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

            #ifdef PARALLEL
            std::copy(POLICY, vector.elements.begin(), vector.elements.end(), this->elements.begin());
            #else
            std::copy(vector.elements.begin(), vector.elements.end(), this->elements.begin());
            #endif

            return *this;
        }

        /**
         * @brief Scalar copy operator.
         * 
         * @param scalar 
         * @return Vector& 
         */
        Vector &operator =(const T &scalar) {
            #ifdef PARALLEL
            std::for_each(POLICY, this->elements.begin(), this->elements.end(), [scalar](auto &element){ element = scalar; });
            #else
            std::for_each(this->elements.begin(), this->elements.end(), [scalar](auto &element){ element = scalar; });
            #endif

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

        // COMPARISONS.

        /**
         * @brief Vector < Scalar.
         * 
         */
        Mask operator <(const T &scalar) const {
            Mask mask(this->length, false);

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elemen, mask.begin(), [scalar](const auto &element){ return element < scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), mask.begin(), [scalar](const auto &element){ return element < scalar; });
            #endif

            return mask;
        }

        /**
         * @brief Vector > Scalar.
         * 
         */
        Mask operator >(const T &scalar) const {
            Mask mask(this->length, false);

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elemen, mask.begin(), [scalar](const auto &element){ return element > scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), mask.begin(), [scalar](const auto &element){ return element > scalar; });
            #endif

            return mask;
        }

        /**
         * @brief Vector < Vector.
         * 
         * @param vector 
         * @return Mask 
         */
        Mask operator <(const Vector &vector) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->length == vector.length);
            #endif

            Mask mask(this->length, false);

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), mask.begin(), [](const auto &first, const auto &second){ return first < second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), mask.begin(), [](const auto &first, const auto &second){ return first < second; });
            #endif

            return mask;
        }

        /**
         * @brief Vector > Vector.
         * 
         * @param vector 
         * @return Mask 
         */
        Mask operator >(const Vector &vector) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->length == vector.length);
            #endif

            Mask mask(this->length, false);

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), mask.begin(), [](const auto &first, const auto &second){ return first > second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), mask.begin(), [](const auto &first, const auto &second){ return first > second; });
            #endif

            return mask;
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
            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [](const auto &element){ return -element; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [](const auto &element){ return -element; });
            #endif

            return result;
        }

        /**
         * @brief Scalar product.
         * 
         * @param scalar 
         * @return Vector 
         */
        Vector operator *(const T &scalar) const {
            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

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
            Vector result{vector.length};

            #ifdef PARALLEL
            std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

            return result;
        }

        /**
         * @brief Scalar product and assignation.
         * 
         * @param scalar 
         * @return Vector& 
         */
        Vector &operator *=(const T &scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element * scalar; });
            #endif

            return *this;
        }

        /**
         * @brief Scalar division.
         * 
         * @param scalar 
         * @return Vector 
         */
        Vector operator /(const T &scalar) const {
            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #endif

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
            Vector result{vector.length};

            #ifdef PARALLEL
            std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return scalar / element; });
            #else
            std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return scalar / element; });
            #endif

            return result;
        }

        /**
         * @brief Scalar division and assignation.
         * 
         * @param scalar 
         * @return Vector& 
         */
        Vector &operator /=(const T &scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element / scalar; });
            #endif

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

            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #endif

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

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first + second; });
            #endif

            return *this;
        }

        /**
         * @brief Scalar "sum".
         * 
         * @param scalar 
         * @return Vector 
         */
        Vector operator +(const T &scalar) const {
            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #endif
            
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
            Vector result{vector.length};

            #ifdef PARALLEL
            std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #else
            std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #endif
            
            return result;
        }

        /**
         * @brief Scalar "sum" and assignation.
         * 
         * @param scalar 
         * @return Vector& 
         */
        Vector &operator +=(const T &scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element + scalar; });
            #endif
            
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

            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #endif

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

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first - second; });
            #endif

            return *this;
        }

        /**
         * @brief Scalar "difference".
         * 
         * @param scalar 
         * @return Vector 
         */
        Vector operator -(const T &scalar) const {
            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element - scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), result.elements.begin(), [scalar](const auto &element){ return element - scalar; });
            #endif
            
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
            Vector result{vector.length};

            #ifdef PARALLEL
            std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return scalar - element; });
            #else
            std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [scalar](const auto &element){ return scalar - element; });
            #endif
            
            return result;
        }

        /**
         * @brief Scalar "difference" and assignation.
         * 
         * @param scalar 
         * @return Vector& 
         */
        Vector &operator -=(const T &scalar) {
            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element - scalar; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), this->elements.begin(), [scalar](const auto &element){ return element - scalar; });
            #endif
            
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

            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first * second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first * second; });
            #endif

            return result;
        }

        /**
         * @brief Vector element-wise product and assignation.
         * 
         * @param vector 
         * @return T 
         */
        Vector operator *=(const Vector &vector) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->length == vector.length);
            #endif

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first * second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first * second; });
            #endif

            return *this;
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

            Vector result{this->length};

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first / second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), result.elements.begin(), [](const auto &first, const auto &second){ return first / second; });
            #endif

            return result;
        }

        /**
         * @brief Vector element-wise division and assignation.
         * 
         * @param vector 
         * @return T 
         */
        Vector operator /=(const Vector &vector) const {
            #ifndef NDEBUG // Integrity check.
            assert(this->length == vector.length);
            #endif

            #ifdef PARALLEL
            std::transform(POLICY, this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first / second; });
            #else
            std::transform(this->elements.begin(), this->elements.end(), vector.elements.begin(), this->elements.begin(), [](const auto &first, const auto &second){ return first / second; });
            #endif

            return *this;
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