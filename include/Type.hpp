/**
 * @file Type.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TYPE_PACS
#define TYPE_PACS

// Concepts.
#include <concepts>

// Math.
#include <complex>
#include <cmath>

// Zero tolerance.
#ifndef TOLERANCE_PACS
#define TOLERANCE_PACS 1E-10
#endif

namespace pacs {

    // Numeric type concept for Matrix and Vector.

    /**
     * @brief Addable types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Addable = requires(T first, T second) {
        {first += second} -> std::convertible_to<T>;
        {first -= second} -> std::convertible_to<T>;
    };

    /**
     * @brief Multipliable types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Multipliable = requires(T first, T second) {
        {first *= second} -> std::convertible_to<T>;
        {first /= second} -> std::convertible_to<T>;
    };

    /**
     * @brief Absolute value supported types.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Absolute = requires(T value) {
        {std::abs(value)} -> std::convertible_to<double>;
    };

    /**
     * @brief Numeric type.
     * 
     * @tparam T 
     */
    template<typename T>
    concept NumericType = Addable<T> && Multipliable<T> && Absolute<T>;

}

#endif