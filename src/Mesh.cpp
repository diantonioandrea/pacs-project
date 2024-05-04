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
    std::vector<Polygon> voronoi(const Polygon &domain, const std::size_t &cells) {
        std::vector<Point> points;

        double x_min = 0.0, x_max = 0.0;
        double y_min = 0.0, y_max = 0.0;

        for(const auto &point: domain.vertices()) {
            if(x_min > point[0])
                x_min = point[0];
            
            if(x_max < point[0])
                x_max = point[0];
            
            if(y_min > point[1])
                y_min = point[1];
            
            if(y_max < point[1])
                y_max = point[1];
        }

        // Padding.
        double x_padding = (x_max - x_min) / 100.0;
        double y_padding = (y_max - y_min) / 100.0;

        x_min += x_padding;
        x_max -= x_padding;

        y_min += y_padding;
        y_max -= y_padding;

        // Points generation.
        for(std::size_t j = 0; j < cells; ++j) {
            double x, y;

            do {
                x = x_min + (x_max - x_min) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
                y = y_min + (y_max - y_min) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
            } while(!(domain.contains({x, y})));

            points.emplace_back(x, y);
        }

        return voronoi(domain, points);
    }

    // LLOYD.

    /**
     * @brief Relaxes a Voronoi diagram through Lloyd's algorithm.
     * 
     * @param domain 
     * @param cells 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> lloyd(const Polygon &domain, const std::vector<Polygon> &cells, const std::size_t &iterations) {
        std::vector<Point> centroids;

        for(const auto &cell: cells)
            centroids.emplace_back(cell.centroid());
        
        if(iterations)
            return lloyd(domain, voronoi(domain, centroids), iterations - 1);

        return voronoi(domain, centroids);
    }

}