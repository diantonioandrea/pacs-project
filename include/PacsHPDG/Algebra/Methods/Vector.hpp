/**
 * @file Vector.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Vector operations and utilities.
 * @date 2024-05-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef VECTOR_METHODS_PACS
#define VECTOR_METHODS_PACS

#include "../Vector.hpp"

namespace pacs {

    /**
     * @brief Get the minimum value from a vector.
     * 
     * @param vector Input vector.
     * @return T 
     */
    template<NumericType T>
    inline T min(const Vector<T> &vector) {
        #ifdef PARALLEL
        return *std::min_element(POLICY, vector.elements.begin(), vector.elements.end());
        #else
        return *std::min_element(vector.elements.begin(), vector.elements.end());
        #endif
    }

    /**
     * @brief Get the maximum value from a vector.
     * 
     * @param vector Input vector.
     * @return T 
     */
    template<NumericType T>
    inline T max(const Vector<T> &vector) {
        #ifdef PARALLEL
        return *std::max_element(POLICY, vector.elements.begin(), vector.elements.end());
        #else
        return *std::max_element(vector.elements.begin(), vector.elements.end());
        #endif
    }

    /**
     * @brief Calculate the sum of all elements in a vector.
     * 
     * @param vector Input vector.
     * @return T 
     */
    template<NumericType T>
    inline T sum(const Vector<T> &vector) {
        #ifdef PARALLEL
        return std::reduce(POLICY, vector.elements.begin(), vector.elements.end(), static_cast<T>(0));
        #else
        return std::reduce(vector.elements.begin(), vector.elements.end(), static_cast<T>(0));
        #endif
    }

    /**
     * @brief Compute the dot product of two vectors.
     * 
     * @param first First vector.
     * @param second Second vector.
     * @return T 
     */
    template<NumericType T>
    inline T dot(const Vector<T> &first, const Vector<T> &second) {
        #ifndef NDEBUG
        assert(first.length == second.length); // Integrity check.
        #endif

        if constexpr (Conjugable<T>) {
            #ifdef PARALLEL
            return std::transform_reduce(POLICY, first.elements.begin(), first.elements.end(), second.elements.begin(), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * std::conj(second); });
            #else
            return std::transform_reduce(first.elements.begin(), first.elements.end(), second.elements.begin(), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * std::conj(second); });
            #endif
        }

        #ifdef PARALLEL
        return std::transform_reduce(POLICY, first.elements.begin(), first.elements.end(), second.elements.begin(), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * second; });
        #else
        return std::transform_reduce(first.elements.begin(), first.elements.end(), second.elements.begin(), static_cast<T>(0), std::plus{}, [](const auto &first, const auto &second){ return first * second; });
        #endif
    }

    /**
     * @brief Compute the L2 norm of a vector.
     * 
     * @param vector Input vector.
     * @return Real 
     */
    template<NumericType T>
    inline Real norm(const Vector<T> &vector) {
        #ifdef PARALLEL
        return std::sqrt(std::transform_reduce(POLICY, vector.elements.begin(), vector.elements.end(), static_cast<T>(0), std::plus{}, [](const auto &element){ return std::abs(element) * std::abs(element); }));
        #else
        return std::sqrt(std::transform_reduce(vector.elements.begin(), vector.elements.end(), static_cast<T>(0), std::plus{}, [](const auto &element){ return std::abs(element) * std::abs(element); }));
        #endif
    }

    /**
     * @brief Concatenate two vectors.
     * 
     * @param first First vector.
     * @param second Second vector.
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> stack(const Vector<T> &first, const Vector<T> &second) {
        Vector<T> result{first.length + second.length};

        for(std::size_t j = 0; j < first.length; ++j)
            result[j] = first[j];

        for(std::size_t j = first.length; j < first.length + second.length; ++j)
            result[j] = second[j - first.length];

        return result;
    }

    /**
     * @brief Reverse a vector.
     * 
     * @param vector Input vector.
     * @return Vector<T> 
     */
    template<NumericType T>
    Vector<T> flip(const Vector<T> &vector) {
        Vector<T> result{vector.length};

        for(std::size_t j = 0; j < result.length; ++j)
            result[j] = vector[vector.length - j - 1];

        return result;
    }

    /**
     * @brief Create a mask for the highest elements in a vector.
     * 
     * @param vector Input vector.
     * @param number Number of elements to mask.
     * @return Mask 
     */
    template<NumericType T>
    Mask highest(const Vector<T> &vector, const std::size_t &number) {
        #ifndef NDEBUG // Integrity check.
        assert(number <= vector.length);
        #endif

        Mask mask(vector.length, false);

        std::vector<T> elements;
        std::vector<std::size_t> indices;

        for(std::size_t j = 0; j < vector.length; ++j) {
            elements.emplace_back(vector[j]);
            indices.emplace_back(j);
        }

        for(std::size_t j = 0; j < elements.size(); ++j) 
            for(std::size_t k = 0; k < elements.size() - 1; ++k)
                if(elements[k] > elements[k + 1]) {
                    std::swap(elements[k], elements[k + 1]);
                    std::swap(indices[k], indices[k + 1]);
                }

        for(std::size_t j = vector.length - number; j < vector.length; ++j)
            mask[indices[j]] = true;

        return mask;
    }

    /**
     * @brief Create a mask for the lowest elements in a vector.
     * 
     * @param vector Input vector.
     * @param number Number of elements to mask.
     * @return Mask 
     */
    template<NumericType T>
    Mask lowest(const Vector<T> &vector, const std::size_t &number) {
        #ifndef NDEBUG // Integrity check.
        assert(number <= vector.length);
        #endif
        
        Mask mask(vector.length, false);

        std::vector<T> elements;
        std::vector<std::size_t> indices;

        for(std::size_t j = 0; j < vector.length; ++j) {
            elements.emplace_back(vector[j]);
            indices.emplace_back(j);
        }

        for(std::size_t j = 0; j < elements.size(); ++j) 
            for(std::size_t k = 0; k < elements.size() - 1; ++k)
                if(elements[k] > elements[k + 1]) {
                    std::swap(elements[k], elements[k + 1]);
                    std::swap(indices[k], indices[k + 1]);
                }

        for(std::size_t j = 0; j < number; ++j)
            mask[indices[j]] = true;

        return mask;
    }
}

namespace std {

    /**
     * @brief Overload of std::abs for vectors.
     * 
     * @param vector Input vector.
     * @return pacs::Vector<T> 
     */
    template<pacs::NumericType T>
    pacs::Vector<T> abs(const pacs::Vector<T> vector) {
        pacs::Vector<T> result{vector.length};

        #ifdef PARALLEL
        std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [](auto &element){ return std::abs(element); });
        #else
        std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [](auto &element){ return std::abs(element); });
        #endif

        return result;
    }

}

#endif