/**
 * @file polygon.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Testing Polygon.
#include <Geometry.hpp>

int main() {

    // Constructing some points.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};
    pacs::Point e{0.5, 0.5};
    pacs::Point f{2.5, 2.5};

    // Constructing a Polygon.
    pacs::Polygon polygon{{a, b, c, d}};

    // Polygon output.
    std::cout << polygon << std::endl;

    // Segments output.
    for(const auto &segment: polygon.edges())
        std::cout << segment << std::endl;
    
    // Point checks.
    std::cout << polygon.contains(e) << std::endl;
    std::cout << polygon.contains(f) << std::endl;

    // Reduction.
    pacs::Line reduction{-1.0, 1.0, -0.5};
    std::cout << pacs::reduce(polygon, reduction, e) << std::endl;

    // Collapse.
    std::cout << pacs::collapse(polygon, a) << std::endl; // Vertex collapse.
    std::cout << pacs::collapse(polygon, {a, b}) << std::endl; // Edge collapse.
    
}