/**
 * @file lloyd.cpp
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
#include <Voronoi.hpp>

int main() {

    // Constructs a domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{1.0, -1.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    
    // Voronoi diagram.
    std::vector<pacs::Polygon> diagram = pacs::voronoi(domain, 20);

    // Relaxation.
    for(std::size_t j = 0; j < 15; ++j)
        diagram = pacs::lloyd(domain, diagram);

    for(const auto &cell: diagram)
        std::cout << cell << std::endl;

}