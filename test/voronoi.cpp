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

    // Some points.
    pacs::Point e{0.5, 0.5};
    pacs::Point f{-0.4, 0.25};
    pacs::Point g{0.0, -0.65};
    pacs::Point o{0.0, 0.0};
    
    // Voronoi diagram.
    std::vector<pacs::Polygon> diagram = pacs::voronoi(domain, {e, f, g});

    for(const auto &cell: diagram)
        std::cout << cell << std::endl;

    // pacs::Line line{-1.0, 1.0, -1.0};
    // std::cout << pacs::reduce(domain, line, o) << std::endl;
}