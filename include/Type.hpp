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
#ifndef TOLERANCE
#define TOLERANCE 1E-14
#endif

// STL Parallelism.
#ifndef _LIBCPP_HAS_NO_INCOMPLETE_PSTL
#define PARALLEL
#include<execution>
using POLICY = std::execution::par_unseq;
#endif

// OpenMP Parallelism.
#ifdef _OPENMP
#include <omp.h>
#endif

namespace pacs {

    // Real alias.

    using Real = long double;

    // Function alias.
    
    using Function = Real (*) (const Real &, const Real &);

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
        {std::abs(value)} -> std::convertible_to<Real>;
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