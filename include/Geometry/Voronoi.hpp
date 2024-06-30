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

#include <Type.hpp>

#include "Shapes.hpp"

#include <vector>

#ifndef GEOMETRY_SAFE
#define GEOMETRY_SAFE 1E-6
#endif

#ifndef GEOMETRY_PADDING
#define GEOMETRY_PADDING 5E-2
#endif

namespace pacs {

    // VORONOI.

    std::vector<Polygon> voronoi(const Polygon &, const std::vector<Point> &, const bool &reflect = false);
    std::vector<Polygon> voronoi_random(const Polygon &, const std::size_t &, const bool &reflect = false);
    std::vector<Polygon> voronoi_uniform(const Polygon &, const std::size_t &, const bool &reflect = false);

    // TRIANGLES.

    std::vector<Polygon> triangulate(const Polygon &);
    std::vector<Polygon> triangulate(const std::vector<Polygon> &);

}

#endif