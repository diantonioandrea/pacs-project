/**
 * @file lshape_domain.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

#include <cmath>
#include <string>

int main() {

    // Domain.
    pacs::Point a{-1.0, -1.0};
    pacs::Point b{0.0, -1.0};
    pacs::Point c{0.0, 0.0};
    pacs::Point d{1.0, 0.0};
    pacs::Point e{1.0, 1.0};
    pacs::Point f{-1.0, 1.0};

    pacs::Polygon domain{{a, b, c, d, e, f}};

    // Meshes.
    for(std::size_t j = 0; j < 4; ++j) {
        
        // Elements.
        std::size_t elements = 100 * std::pow(2, j);

        // Mesh.
        pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, elements, false, true)};

        // Output.
        mesh.write("output/lshape_" + std::to_string(elements) + ".poly");
    }

}