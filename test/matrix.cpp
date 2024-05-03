/**
 * @file matrix.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include<iostream>

// Testing Matrix.
#include <Matrix.hpp>

int main() {

    // Constructing a matrix.
    pacs::Matrix<double> matrix{2, 2};

    // Write.
    matrix(0, 0) = 1;
    matrix(1, 1) = -1;
    
    // Output.
    std::cout << matrix << std::endl;

    // Vector product.
    pacs::Vector<double> vector{2};
    
    vector[0] = 1;
    vector[1] = 2;

    // Vector product output.
    std::cout << (matrix * vector) << std::endl;
}