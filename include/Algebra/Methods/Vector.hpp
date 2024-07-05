/**
 * @file Vector.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
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
     * @brief Minimum value of a Vector.
     * 
     * @param vector 
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
     * @brief Maximum value of a Vector.
     * 
     * @param vector 
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
     * @brief Sum of the elements inside a Vector.
     * 
     * @tparam T 
     * @param vector 
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
     * @brief Vector dot product.
     * 
     * @param first 
     * @param second 
     * @return T 
     */
    template<NumericType T>
    inline T dot(const Vector<T> &first, const Vector<T> &second) {
        #ifndef NDEBUG
        assert(first.length == second.length);
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
     * @brief Returns the l2 norm of the Vector.
     * 
     * @param vector 
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

        for(std::size_t j = 0; j < first.length; ++j)
            result[j] = first[j];

        for(std::size_t j = first.length; j < first.length + second.length; ++j)
            result[j] = second[j - first.length];

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
    
    /**
     * @brief Marks the highest elements in a Vector.
     * 
     * @tparam T 
     * @param vector 
     * @param number 
     * @return Mask 
     */
    template<NumericType T>
    Mask highest(const Vector<T> &vector, const std::size_t &number) {
        #ifndef NDEBUG // Integrity check.
        assert(number <= vector.length);
        #endif

        Mask mask(vector.length, false);

        // Indices and copy of vector.
        std::vector<T> elements;
        std::vector<std::size_t> indices;

        for(std::size_t j = 0; j < vector.length; ++j) {
            elements.emplace_back(vector[j]);
            indices.emplace_back(j);
        }

        // Sorting.
        for(std::size_t j = 0; j < elements.size(); ++j) 
            for(std::size_t k = 0; k < elements.size() - 1; ++k)
                if(elements[k] > elements[k + 1]) {
                    std::swap(elements[k], elements[k + 1]);
                    std::swap(indices[k], indices[k + 1]);
                }

        // Masking.
        for(std::size_t j = vector.length - number; j < vector.length; ++j)
            mask[indices[j]] = true;

        return mask;
    }

    /**
     * @brief Marks the lowest elements in a Vector.
     * 
     * @tparam T 
     * @param vector 
     * @param number 
     * @return Mask 
     */
    template<NumericType T>
    Mask lowest(const Vector<T> &vector, const std::size_t &number) {
        #ifndef NDEBUG // Integrity check.
        assert(number <= vector.length);
        #endif
        
        Mask mask(vector.length, false);

        // Indices and copy of vector.
        std::vector<T> elements;
        std::vector<std::size_t> indices;

        for(std::size_t j = 0; j < vector.length; ++j) {
            elements.emplace_back(vector[j]);
            indices.emplace_back(j);
        }

        // Sorting.
        for(std::size_t j = 0; j < elements.size(); ++j) 
            for(std::size_t k = 0; k < elements.size() - 1; ++k)
                if(elements[k] > elements[k + 1]) {
                    std::swap(elements[k], elements[k + 1]);
                    std::swap(indices[k], indices[k + 1]);
                }

        // Masking.
        for(std::size_t j = 0; j < number; ++j)
            mask[indices[j]] = true;

        return mask;
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

        #ifdef PARALLEL
        std::transform(POLICY, vector.elements.begin(), vector.elements.end(), result.elements.begin(), [](auto &element){ return std::abs(element); });
        #else
        std::transform(vector.elements.begin(), vector.elements.end(), result.elements.begin(), [](auto &element){ return std::abs(element); });
        #endif

        return result;
    }

}

#endif