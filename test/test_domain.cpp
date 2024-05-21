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
#include <iostream>

// Containers.
#include <vector>

// Testing Mesh.
#include <Geometry.hpp>
#include <Mesh.hpp>

int main() {

    // Constructs a domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};
    
    // Constructing a mesh.
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, 64, true)};

    // Mesh output.
    mesh.write("output/mesh.poly");

}