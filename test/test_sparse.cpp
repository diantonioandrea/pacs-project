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
#include <iostream>

// Testing Sparse.
#include <Sparse.hpp>
using pacs::Real;

int main() {

    // Constructing a matrix.
    pacs::Sparse<Real> sparse{2, 2};
    
    // Insert.
    sparse.insert(0, 0, 1);
    sparse.insert(1, 1, -1);

    // Compression.
    sparse.compress();

    // Output.
    std::cout << sparse << std::endl;

    // Vector product.
    pacs::Vector<Real> vector{2};
    
    vector[0] = 1;
    vector[1] = 2;

    // Vector product output.
    std::cout << (sparse * vector) << std::endl;
    
}