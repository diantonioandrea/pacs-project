/**
 * @file voronoi.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include<iostream>

// Containers.
#include <vector>

// Testing Voronoi.
#include <Geometry.hpp>
#include <Mesh.hpp>

int main() {

    // Constructs a domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{1.0, -1.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Picks some points.
    pacs::Point e{0.5, 0.0};
    pacs::Point f{-0.5, 0.5};
    pacs::Point g{0.0, -0.75};

    // Cells output.
    for(const auto &cell: pacs::voronoi(domain, {e, f, g}))
        std::cout << cell << std::endl;
    
}