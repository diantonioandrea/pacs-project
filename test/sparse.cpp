/**
 * @file sparse.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include<iostream>

// Testing Sparse.
#include <Sparse.hpp>

int main() {

    // Constructing a matrix.
    pacs::Sparse<double> sparse{2, 2};
    
    // Insert.
    sparse.insert(0, 0, 1.0);
    sparse.insert(1, 1, 1.0);

    // Compression.
    sparse.compress();

    // Output.
    std::cout << sparse << std::endl;

    // Vector product.
    pacs::Vector<double> vector{2};
    
    vector[0] = 1;
    vector[1] = -1;

    // Vector product output.
    std::cout << (sparse * vector) << std::endl;
}