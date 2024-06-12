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
    pacs::Mesh mesh{domain, pacs::mesh_diagram("data/square_30.poly")};

    // Refining the first four elements.
    pacs::Mask refinement(30);

    refinement[0] = true;
    refinement[1] = true;
    refinement[2] = true;

    // Constructing a new mesh.
    std::vector<pacs::Polygon> refined = pacs::mesh_refine_size(mesh, refinement);
    pacs::Mesh refined_mesh{domain, refined};

    // Mesh output.
    refined_mesh.write("output/refined.poly");

}