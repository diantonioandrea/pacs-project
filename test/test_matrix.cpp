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
#include <iostream>

// Testing Matrix.
#include <Matrix.hpp>
using pacs::Real;

int main() {

    // Constructing a matrix.
    pacs::Matrix<Real> matrix{2, 2};

    // Write.
    matrix(0, 0) = 1;
    matrix(1, 1) = -1;
    
    // Output.
    std::cout << matrix << std::endl;

    // Vector product.
    pacs::Vector<Real> vector{2};
    
    vector[0] = 1;
    vector[1] = 2;

    // Vector product output.
    std::cout << (matrix * vector) << std::endl;

    // .row() and .column().
    std::cout << matrix.row(0) << std::endl;
    matrix.column(0, 3.0);
    std::cout << matrix << std::endl << std::endl;

    // Product output.
    std::cout << matrix * matrix << std::endl << std::endl;

    // Transpose output.
    std::cout << matrix.transpose() << std::endl;
    
}