/**
 * @file Voronoi.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef VORONOI_PACS
#define VORONOI_PACS

// Type.
#include <Type.hpp>

// Geometry.
#include <Geometry.hpp>

// Containers.
#include <vector>

namespace pacs {

    // VORONOI.

    std::vector<Polygon> voronoi(const Polygon &, const std::vector<Point> &, const bool &reflect = false);
    std::vector<Polygon> voronoi(const Polygon &, const std::size_t &, const bool &reflect = false);

    // TRIANGLES.

    std::vector<Polygon> triangulate(const Polygon &);
    std::vector<Polygon> triangulate(const std::vector<Polygon> &);

}

#endif