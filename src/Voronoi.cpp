/**
 * @file Voronoi.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Voronoi.hpp>

// Random.
#include <cstdlib>
#include <ctime>

// Safe distance.
#ifndef GEOMETRY_SAFE
#define GEOMETRY_SAFE 0.05
#endif

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

                // Failure.
                if(!cell.contains(point))
                    return std::vector<Polygon>{};

                Line line = bisector(point, points[k]);
                cell = reduce(cell, line, point);
            }

            cells.emplace_back(cell);
        }

        return cells;
    }

    /**
     * @brief Returns the Voronoi diagram of a random number of points inside a bounded domain.
     * 
     * @param domain 
     * @param cells 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi(const Polygon &domain, const std::size_t &cells) {
        // Seeding.
        std::srand(time(0));

        std::vector<Point> points;
        std::size_t counter = 0;

        do {
            Point random = domain.random();
            bool check = true;

            for(const auto &point: points)
                if((point - random).norm() <= GEOMETRY_SAFE) {
                    check = false;
                    break;
                }

            if(check) {
                points.emplace_back(random);
                ++counter;
            }

        } while(counter < cells);

        std::vector<Polygon> diagram = voronoi(domain, points);

        if(diagram.size())
            return diagram;

        // Retries.
        return voronoi(domain, cells);
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

        std::vector<Polygon> diagram = voronoi(domain, centroids);

        if(diagram.size())
            return diagram;

        // Failure.
        return cells;
    }

    // TRIANGLES.

    std::vector<Polygon> triangulate(const Polygon &polygon) {
        std::vector<Polygon> triangles;
        Point centroid = polygon.centroid();

        for(const auto &edge: polygon.edges()) {
            std::vector<Point> vertices{edge[0], edge[1], centroid};
            triangles.emplace_back(vertices);
        }

        return triangles;
    }

    std::vector<Polygon> triangulate(const std::vector<Polygon> &polygons) {
        std::vector<Polygon> triangles;

        for(const auto &polygon: polygons) {
            for(const auto &triangle: triangulate(polygon))
                triangles.emplace_back(triangle);
        }

        return triangles;
    }

}