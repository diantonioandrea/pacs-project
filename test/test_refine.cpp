/**
 * @file test_refine.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

#include <iostream>
#include <vector>

int main() {

    // Constructs a domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{1.0, -1.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};
    
    // Constructing a mesh.
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, 32)};

    // Refining the first four elements.
    std::vector<std::size_t> indices{0, 1, 2, 3};
    std::vector<pacs::Polygon> refined = pacs::mesh_refine(mesh, indices);

    // Constructing a new mesh.
    pacs::Mesh refined_mesh{domain, refined};

    // Mesh output.
    refined_mesh.write("output/refined.poly");

}