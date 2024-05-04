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
#include<iostream>

// Testing Vector.
#include <Vector.hpp>

int main() {

    // Constructing a vector.
    pacs::Vector<double> vector{4};

    // Output.
    std::cout << vector << std::endl;
    
}