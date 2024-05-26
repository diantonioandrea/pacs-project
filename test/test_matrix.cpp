/**
 * @file matrix.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Algebra.hpp>

#include <iostream>

int main() {

    // Constructing a matrix.
    pacs::Matrix<pacs::Real> matrix{3, 2};

    // Write.
    matrix(0, 0) = 1;
    matrix(0, 1) = 2;
    matrix(1, 0) = 3;
    matrix(1, 1) = 4;
    matrix(2, 0) = 5;
    matrix(2, 1) = 6;
    
    // QR decomposition.
    auto [Q, R] = QR(matrix.transpose());

    // Output.
    std::cout << matrix.transpose() << std::endl << std::endl;
    std::cout << Q << std::endl << std::endl;
    std::cout << R << std::endl << std::endl;
    std::cout << Q * Q.transpose() << std::endl << std::endl;
    std::cout << Q * R << std::endl;
}