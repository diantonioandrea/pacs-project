// clang-format off
/**
 * @file Voronoi.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <ranges>

namespace pacs {

    // VORONOI.

    /**
     * @brief Returns the Voronoi diagram of a vector of points inside a bounded domain.
     * 
     * @param domain Polygonal domain.
     * @param points Points.
     * @param reflect Reflection flag.
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi(const Polygon &domain, const std::vector<Point> &points, const bool &reflect) {
        std::vector<Polygon> cells;
        Polygon box = domain;

        // Reflections.
        std::vector<std::vector<Point>> reflection_points;

        if(reflect) {
            for(const auto &point: points) {
                std::vector<Point> reflection_set;

                for(const auto &reflection: reflections(domain, point))
                    reflection_set.emplace_back(reflection);

                reflection_points.emplace_back(reflection_set);
            }

            auto [xy_min, xy_max] = domain.box();
            Real x_min = xy_min[0], y_min = xy_min[1];
            Real x_max = xy_max[0], y_max = xy_max[1];
            
            // std library provides std::minmax_element
            for(const auto &reflection: reflection_points)
                for(const auto &point: reflection) {
                    x_min = (point[0] < x_min) ? point[0] : x_min;
                    x_max = (point[0] > x_max) ? point[0] : x_max;
                    y_min = (point[1] < y_min) ? point[1] : y_min;
                    y_max = (point[1] > y_max) ? point[1] : y_max;
                }

            x_min -= GEOMETRY_PADDING;
            x_max += GEOMETRY_PADDING;
            y_min -= GEOMETRY_PADDING;
            y_max += GEOMETRY_PADDING;

            // Boundary rectangle {a, b, c, d};
            Point a{x_min, y_min};
            Point b{x_max, y_min};
            Point c{x_max, y_max};
            Point d{x_min, y_max};

            box = Polygon{{a, b, c, d}};
        }

        // Cells initialization.
        cells.resize(points.size(), box);

        #pragma omp parallel for
        for(std::size_t j = 0; j < points.size(); ++j) {
            for(std::size_t k = 0; k < points.size(); ++k) {
                if(k == j)
                    continue;

                cells[j] = reduce(cells[j], bisector(points[j], points[k]), points[j]);
            }

            if(reflect)
                for(const auto &reflection: reflection_points[j])
                    cells[j] = reduce(cells[j], bisector(points[j], reflection), points[j]);
        }

        return cells;
    }

    /**
     * @brief Returns the Voronoi diagram of a random number of points inside a bounded domain.
     * 
     * @param domain Polygonal domain.
     * @param cells Number of points.
     * @param reflect Reflection flag.
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi_random(const Polygon &domain, const std::size_t &cells, const bool &reflect) {

        // Seeding.
        std::srand(time(0));

        std::vector<Point> points;
        std::size_t counter = 0;

        do {
            Point random = domain.random();
            bool check = true;

            for(const auto &point: points) {
                if(distance(point, random) <= GEOMETRY_SAFE) {
                    check = false;
                    break;
                }
            }

            if(check) {
                points.emplace_back(random);
                ++counter;
            }

        } while(counter < cells);

        return voronoi(domain, points, reflect);
    }

    /**
     * @brief Returns the Voronoi diagram of uniformely distributes points inside a bounded domain.
     * 
     * @param domain Polygonal domain.
     * @param cells Number of points.
     * @param reflect Reflection flag.
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> voronoi_uniform(const Polygon &domain, const std::size_t &cells, const bool &reflect) {

        // Points.
        std::vector<Point> points;

        // Grid.
        std::size_t cells_x = std::sqrt(cells);
        std::size_t cells_y = cells / cells_x;

        // Bounds.
        auto [xy_min, xy_max] = domain.box();
        Real y_min = xy_min[1] + GEOMETRY_PADDING, y_max = xy_max[1] - GEOMETRY_PADDING;

        // Y step.
        Real y_step = (y_max - y_min) / static_cast<Real>(cells_y);

        for(std::size_t j = 0; j < cells_y; ++j) {
            std::vector<Point> bounds = intersections(Line{0.0, 1.0, y_min + j * y_step}, domain);

            // Start and end.
            Point start = bounds[0];
            Point end = bounds[0];

            for(const auto &bound: bounds) {
                start = (bound[0] < start[0]) ? bound : start;
                end = (bound[0] > end[0]) ? bound : end;
            }

            start[0] += GEOMETRY_PADDING;
            end[0] -= GEOMETRY_PADDING;

            // X step.
            Real x_step = (end[0] - start[0]) / static_cast<Real>(cells_x);

            // Points.
            for(std::size_t k = 0; k < cells_x; ++k) {
                Point point{start[0] + k * x_step, y_min + j * y_step};

                if(domain.contains(point))
                    points.emplace_back(point);
            }
        }

        return voronoi(domain, points, reflect);
    }

    // TRIANGLES.

    /**
     * @brief Triangulates a polygon.
     * 
     * @param polygon Polygon.
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> triangulate(const Polygon &polygon) {
        std::vector<Polygon> triangles;
        Point centroid = polygon.centroid();

        for(const auto &edge: polygon.edges()) {
            std::vector<Point> vertices{edge[0], edge[1], centroid};
            triangles.emplace_back(vertices);
        }

        return triangles;
    }
    
    /**
     * @brief Triangulates many polygons.
     * 
     * @param polygons Polygons.
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> triangulate(const std::vector<Polygon> &polygons) {
        std::vector<Polygon> triangles;

        for(const auto &polygon: polygons) {
            for(const auto &triangle: triangulate(polygon))
                triangles.emplace_back(triangle);
        }

        return triangles;
    }

}