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
using pacs::Real;


int main() {

    // Constructs a Matrix.
    pacs::Matrix<Real> matrix{2, 2};

    matrix(0, 0) = 4.0;
    matrix(0, 1) = 1.0;
    matrix(1, 0) = 1.0;
    matrix(1, 1) = 3.0;

    // Constructs a Vector.
    pacs::Vector<Real> vector{2};

    vector[0] = 1.0;
    vector[1] = 2.0;

    // Linear system (Ax = b) solution.
    std::cout << pacs::solve(matrix, vector) << std::endl;
    
}