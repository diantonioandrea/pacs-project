/**
 * @file triangles.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include <iostream>

// Containers.
#include <vector>

// Testing Polygon triangulation.
#include <Geometry.hpp>
#include <Voronoi.hpp>
using pacs::Real;

int main() {

    // Constructs a domain.
    pacs::Point a{-1.5, -1.0};
    pacs::Point b{1.0, -2.0};
    pacs::Point c{2.0, 0.0};
    pacs::Point d{1.0, 1.0};
    pacs::Point e{-1.0, 1.0};

    pacs::Polygon polygon{{a, b, c, d, e}};

    // Centroid triangulation.
    std::vector<pacs::Polygon> triangles = pacs::triangulate(polygon);

    // Triangulation output.
    for(const auto &triangle: triangles)
        std::cout << triangle << std::endl;
}