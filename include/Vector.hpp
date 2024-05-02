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

            Vector(const std::size_t &length, const std::vector<T> &elements): length{length}, elements{elements} {
                #ifndef NDEBUG // Integrity check.
                assert(length > 0);
                assert(elements.size() == length);
                #endif
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