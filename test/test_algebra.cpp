/**
 * @file algebra.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Testing Algebra.
#include <Algebra.hpp>

int main() {

    // Constructs a Matrix.
    pacs::Sparse<pacs::Real> matrix{2, 2};

    matrix.insert(0, 0, 4.0);
    matrix.insert(0, 1, 1.0);
    matrix.insert(1, 0, 1.0);
    matrix.insert(1, 1, 3.0);

    // Constructs a Vector.
    pacs::Vector<pacs::Real> vector{2};

    vector[0] = 1.0;
    vector[1] = 2.0;

    // Linear system (Ax = b) solution.
    std::cout << pacs::solve(matrix, vector) << std::endl;
    
}