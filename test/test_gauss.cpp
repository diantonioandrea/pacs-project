/**
 * @file gauss.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Testing Quadrature.
#include <Quadrature.hpp>
using pacs::Real;

int main() {

    // Some Gauss-Legendre nodes and weights over [0, 1].
    auto [nodes_x, nodes_y, weights] = pacs::quadrature_2d(3);

    // Output
    std::cout << nodes_x << std::endl; // Nodes_x.
    std::cout << nodes_y << std::endl; // Nodes_y.
    std::cout << weights << std::endl; // Weights.
    
}