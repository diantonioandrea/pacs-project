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
#include <Algebra.hpp>

int main() {

    // Sparse.

    // Constructs a Matrix.
    pacs::Sparse<pacs::Real> sparse{2, 2};

    sparse.insert(0, 0, 4.0);
    sparse.insert(0, 1, 1.0);
    sparse.insert(1, 0, 1.0);
    sparse.insert(1, 1, 3.0);

    sparse.compress();

    // Constructs a Vector.
    pacs::Vector<pacs::Real> vector{2};

    vector[0] = 1.0;
    vector[1] = 2.0;

    // Linear system (Ax = b) solution.
    std::cout << pacs::solve(sparse, vector, pacs::GMRES) << std::endl; // GMRES.

    // Dense.
    
    // Constructs a Matrix.
    pacs::Matrix<pacs::Real> dense{2, 2};

    dense(0, 0) = 4.0;
    dense(0, 1) = 1.0;
    dense(1, 0) = 1.0;
    dense(1, 1) = 3.0;

    // Linear system (Ax = b) solution.
    std::cout << pacs::solve(dense, vector, pacs::QRD) << std::endl; // QR.
    std::cout << pacs::solve(dense, vector, pacs::LUD) << std::endl; // LU.
    
}