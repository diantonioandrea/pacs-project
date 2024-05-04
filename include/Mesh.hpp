/**
 * @file Mesh.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH_PACS
#define MESH_PACS

// Geometry.
#include <Geometry.hpp>

// Containers.
#include <vector>

namespace pacs {

    // VORONOI.

    std::vector<Polygon> voronoi(const Polygon &, const std::vector<Point> &);
    std::vector<Polygon> voronoi(const Polygon &, const std::size_t &);

}

#endif