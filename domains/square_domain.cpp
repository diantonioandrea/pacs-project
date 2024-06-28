/**
 * @file square_domain.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Square domain mesh generator.
 * @date 2024-06-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

#include <cmath>
#include <string>

int main(int argc, char **argv) {

    // Domain.
    pacs::Point a{0.0, 0.0};
    pacs::Point b{1.0, 0.0};
    pacs::Point c{1.0, 1.0};
    pacs::Point d{0.0, 1.0};

    pacs::Polygon domain{{a, b, c, d}};

    // Elements.
    if(argc <= 1) {
        std::cout << "Usage: " << argv[0] << " ELEMENTS." << std::endl;
        std::exit(-1);
    }

    std::size_t elements = static_cast<std::size_t>(std::stoi(argv[1]));

    // Mesh.
    pacs::Mesh mesh{domain, pacs::mesh_diagram(domain, elements)};

    // Output.
    mesh.write("output/square_" + std::to_string(elements) + ".poly");

}