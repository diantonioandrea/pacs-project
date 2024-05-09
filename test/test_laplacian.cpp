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
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, 32)};

    // Writes mesh informations to a file.
    mesh.write("laplacian.poly");

    // Builds the laplacian matrix.
    auto [laplacian, dg_laplacian] = pacs::laplacian(mesh);

    // Matrix output.
    std::cout << laplacian << std::endl;
    
}