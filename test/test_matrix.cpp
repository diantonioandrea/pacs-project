/**
 * @file matrix.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>

int main() {

    // Constructing a matrix.
    pacs::Matrix<pacs::Real> matrix{3, 3};

    // Write.
    matrix(0, 0) = 1;
    matrix(0, 1) = 2;
    matrix(0, 2) = 3;
    matrix(1, 0) = 2;
    matrix(1, 1) = 3;
    matrix(1, 2) = 4;
    matrix(2, 0) = -2;
    matrix(2, 1) = 7;
    matrix(2, 2) = 12;
    
    // QR decomposition.
    auto [Q, R] = pacs::QR(matrix);

    // Output.
    std::cout << matrix << std::endl << std::endl;
    std::cout << Q << std::endl << std::endl;
    std::cout << R << std::endl << std::endl;
    std::cout << Q * Q.transpose() << std::endl;
    std::cout << Q * R << std::endl << std::endl;

    // LU decomposition.
    auto [L, U] = pacs::LU(matrix);

    // Output.
    std::cout << matrix << std::endl << std::endl;
    std::cout << L << std::endl << std::endl;
    std::cout << U << std::endl << std::endl;
    std::cout << L * U << std::endl << std::endl;

    // Inversion.
    auto I = pacs::solve(matrix, pacs::identity<pacs::Real>(3));

    // Output.
    std::cout << matrix * I << std::endl << std::endl;
}