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

// Testing Solvers.
#include <Sparse.hpp>

int main() {

    // Sparse.

    // Constructs a Matrix.
    pacs::Sparse<pacs::Real> matrix{2, 2};

    matrix.insert(0, 0, 4.0);
    matrix.insert(0, 1, 1.0);
    matrix.insert(1, 0, 1.0);
    matrix.insert(1, 1, 3.0);

    matrix.compress();

    // Constructs a Vector.
    pacs::Vector<pacs::Real> vector{2};

    vector[0] = 1.0;
    vector[1] = 2.0;

    // Linear system (Ax = b) solution.
    std::cout << matrix.solve<pacs::CG>(vector) << std::endl; // Conjugate gradient.
    std::cout << matrix.solve<pacs::SD>(vector) << std::endl; // Gradient descent.
    std::cout << matrix.solve<pacs::MR>(vector) << std::endl; // Minimal residual.
    std::cout << matrix.solve<pacs::NS>(vector) << std::endl; // Norm Steepest Descent.
    std::cout << matrix.solve<pacs::GS>(vector) << std::endl; // Gauss-Seidel.
    std::cout << matrix.solve<pacs::RFOM>(vector) << std::endl; // Restarted FOM.
    std::cout << matrix.solve<pacs::KM>(vector) << std::endl; // Kaczmarz.
    std::cout << matrix.solve<pacs::RKM>(vector) << std::endl; // Randomized Kaczmarz.

    // Dense.
    
    // Constructs a Matrix.
    pacs::Matrix<pacs::Real> dense{2, 2};

    dense(0, 0) = 4.0;
    dense(0, 1) = 1.0;
    dense(1, 0) = 1.0;
    dense(1, 1) = 3.0;

    // Linear system (Ax = b) solution.
    std::cout << dense.solve(vector) << std::endl; // LU.
    
}