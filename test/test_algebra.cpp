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

// Testing Sparse solvers.
#include <Sparse.hpp>

int main() {

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
    std::cout << matrix.solve(vector) << std::endl; // Conjugate gradient.
    std::cout << matrix.solve<pacs::Descent>(vector) << std::endl; // Gradient descent.
    std::cout << matrix.solve<pacs::Minimal>(vector) << std::endl; // Minimal residual.
    std::cout << matrix.solve<pacs::Gauss>(vector) << std::endl; // Gauss-Seidel.
    std::cout << matrix.solve<pacs::Kaczmarz>(vector) << std::endl; // Kaczmarz.
    std::cout << matrix.solve<pacs::RandomKaczmarz>(vector) << std::endl; // Randomized Kaczmarz.
    
}