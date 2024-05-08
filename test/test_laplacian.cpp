/**
 * @file test_laplacian.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */


// Output.
#include <iostream>

// Testing Laplacian.
#include <Laplacian.hpp>

int main() {

    // Constructs a mesh.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{1.0, -1.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, 32};

    // Writes mesh informations to a file.
    mesh.write("laplacian.poly");

    // Builds the laplacian matrix.
    pacs::Sparse<double> laplacian = pacs::laplacian(mesh);

    // Matrix output.
    std::cout << laplacian << std::endl;
    
}