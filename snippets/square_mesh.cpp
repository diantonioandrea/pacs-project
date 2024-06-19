/**
 * @file square_mesh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-06-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>
using namespace pacs;

int main() {

    // Domain.
    Point a{0.0, 0.0};
    Point b{1.0, 0.0};
    Point c{1.0, 1.0};
    Point d{0.0, 1.0};

    Polygon domain{{a, b, c, d}};

    // Diagram.
    std::vector<Polygon> diagram = mesh_diagram(domain, 100);

    // Mesh.
    Mesh mesh{domain, diagram};

}