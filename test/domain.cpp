/**
 * @file domain.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include<iostream>

// Containers.
#include <vector>

// Testing Mesh.
#include <Geometry.hpp>
#include <Mesh.hpp>

int main() {

    // Constructs a domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{1.0, -1.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    
    // Constructing a mesh.
    pacs::Mesh mesh{domain, 500};

    // Mesh output.
    mesh.write("mesh.poly");

}