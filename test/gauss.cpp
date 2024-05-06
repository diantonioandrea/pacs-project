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
#include<iostream>

// Testing Quadrature.
#include <Quadrature.hpp>

int main() {

    // Some Gauss-Legendre nodes and weights over [0, 1].
    std::vector<pacs::Vector<double>> gauss_legendre = pacs::gauss_legendre(0, 1, 11);

    // Output
    std::cout << gauss_legendre[0] << std::endl;
    std::cout << gauss_legendre[1] << std::endl;
    
}