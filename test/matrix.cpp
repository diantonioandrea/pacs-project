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

    // Output.
    std::cout << matrix << std::endl;
}