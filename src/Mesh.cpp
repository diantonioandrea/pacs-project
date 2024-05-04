/**
 * @file Mesh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Mesh.hpp>

namespace pacs {

    // VORONOI.

    /**
     * @brief Returns the Voronoi diagram of a vector of points inside a bounded domain.
     * 
     * @param domain 
     * @param points 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi(const Polygon &domain, const std::vector<Point> &points) {
        std::vector<Polygon> cells;

        for(std::size_t j = 0; j < points.size(); ++j) {
            Polygon cell = domain;

            for(std::size_t k = 0; k < points.size(); ++k) {
                if(k == j)
                    continue;

                cell = reduce(cell, bisector(points[j], points[k]), points[j]);
            }

            cells.emplace_back(cell);
        }

        return cells;
    }

    /**
     * @brief Returns the Voronoi diagram of a random set of points.
     * 
     * @param domain 
     * @param points 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi(const Polygon &domain, const std::size_t &points) {

        return std::vector<Polygon>{};

    }

}