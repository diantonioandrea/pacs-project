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

// Random.
#include <cstdlib>
#include <ctime>

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
            Point point{points[j]};
            Polygon cell = domain;

            for(std::size_t k = 0; k < points.size(); ++k) {
                if(k == j)
                    continue;

                Line line = bisector(point, points[k]);
                cell = reduce(cell, line, point);
            }

            cells.emplace_back(cell);
        }

        return cells;
    }

    // LLOYD.

    /**
     * @brief Relaxes a Voronoi diagram through a single-step Lloyd's algorithm.
     * 
     * @param domain 
     * @param cells 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> lloyd(const Polygon &domain, const std::vector<Polygon> &cells) {
        std::vector<Point> centroids;

        for(const auto &cell: cells)
            centroids.emplace_back(cell.centroid());

        return voronoi(domain, centroids);
    }

}