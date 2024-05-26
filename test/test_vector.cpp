/**
 * @file vector.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Testing Vector.
#include <Algebra.hpp>

int main() {

    // Constructing two vectors.
    pacs::Vector<pacs::Real> first{4, 1.0};
    pacs::Vector<pacs::Real> second{4};

    second[1] = 2.0;
    second[2] = -2.0;

    // Operations output.
    std::cout << first << std::endl;
    std::cout << second << std::endl;

    std::cout << first + second << std::endl;
    std::cout << first - second << std::endl;

    std::cout << (first += 2.0) << std::endl;
    std::cout << (second -= 2.0) << std::endl;

    std::cout << (first * 3.0) << std::endl;
    std::cout << (second * 3.0) << std::endl;

    std::cout << first << std::endl;
    std::cout << second << std::endl;

    std::cout << first * second << std::endl;

    std::cout << first << std::endl;
    std::cout << second << std::endl;

}