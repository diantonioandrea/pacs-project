/**
 * @file Geometry_Methods.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

// Assertions.
#include <cassert>

// Math.
#include <cmath>

// Swap.
#include <algorithm>

namespace pacs {

    /**
     * @brief Finds the bisector given two Points.
     * 
     * @param p 
     * @param q 
     * @return Line 
     */
    Line bisector(const Point &p, const Point &q) {
        #ifndef NDEBUG // Checks.
        assert(!(p - q).is_zero());
        #endif

        double px = p[0], qx = q[0];
        double py = p[1], qy = q[1];
        double mx = (px + qx) / 2, my = (py + qy) / 2;

        // Evaluation by cases.

        // Px = Qx.
        if(std::abs(px - qx) <= GEOMETRY_TOLERANCE)
            return Line{0.0, 1.0, my};

        // Py = Qy.
        if(std::abs(py - qy) <= GEOMETRY_TOLERANCE)
            return Line{1.0, 0.0, mx};

        // Default.
        double a = qx - px;
        double b = qy - py;
        double c = a * mx + b * my;
        return Line{a, b, c};
    }

    /**
     * @brief Returns the intersection between two Lines.
     * 
     * @param s 
     * @param r 
     * @return std::vector<Point> 
     */
    std::vector<Point> intersections(const Line &s, const Line &r) {
        std::vector<Point> points;

        // Evaluation by cases.

        // No (or infinite) intersections.
        if(s.is_parallel(r))
            return points;

        // S is vertical.
        if(std::abs(s[1]) <= GEOMETRY_TOLERANCE) {
            points.emplace_back(s[2] / s[0], r.y(s[2] / s[0]));
            return points;
        }

        // R is vertical.
        if(std::abs(r[1]) <= GEOMETRY_TOLERANCE) {
            points.emplace_back(r[2] / r[0], s.y(r[2] / r[0]));
            return points;
        }

        // Default case.

        // s: y = ax + c.
        double a = - s[0] / s[1];
        double c = s[2] / s[1];

        // r: y = bx + d.
        double b = -r[0] / r[1];
        double d = r[2] / r[1];

        double x = (d - c) / (a - b);
        double y = a * x + c;

        points.emplace_back(x, y);
        return points;
    }

    /**
     * @brief Returns the intersection between a Line and a Segment.
     * 
     * @param line 
     * @param segment 
     * @return std::vector<Point> 
     */
    std::vector<Point> intersections(const Line &line, const Segment &segment) {
        std::vector<Point> points;

        for(const auto &point: intersections(line, segment.line())) {
            if(segment.contains(point))
                points.emplace_back(point);
        }

        return points;
    }

    /**
     * @brief Returns the intersection(s) between a line and a Polygon.
     * 
     * @param line 
     * @param polygon 
     * @return std::vector<Point> 
     */
    std::vector<Point> intersections(const Line &line, const Polygon &polygon) {
        std::vector<Point> points;

        for(const auto &segment: polygon.edges()) {
            for(const auto &point: intersections(line, segment))
                points.emplace_back(point);
        }

        return points;
    }

    /**
     * @brief Reduce a Polygon cutted by a Line. Picks the part which contains the argument Point.
     * 
     * @param polygon 
     * @param line 
     * @param point 
     * @return Polygon 
     */
    Polygon reduce(const Polygon &polygon, const Line &line, const Point &point) {
        std::vector<Point> points = intersections(line, polygon);
        std::vector<Point> new_vertices;

        #ifndef NDEBUG // Integrity check.
        assert(polygon.contains(point));
        assert(points.size() <= 2); // Convex domain assumption.
        #endif

        if(points.size() <= 1)
            return polygon;

        // Building.
        std::vector<Point> vertices = polygon.vertices();
        std::vector<Segment> edges = polygon.edges();
        std::size_t index = 0;

        if(line > point) { // Point under the Line.

            if((std::abs(line[0]) <= GEOMETRY_TOLERANCE) && (points[0][0] < points[1][0])) {
                std::swap(points[0], points[1]);
            } else if((line.angular() > 0.0) && (points[0][1] < points[1][1])) {
                std::swap(points[0], points[1]);
            } else if((line.angular() < 0.0) && (points[0][1] > points[1][1])) {
                std::swap(points[0], points[1]);
            }

            new_vertices.emplace_back(points[0]);
            new_vertices.emplace_back(points[1]);

            for(std::size_t j = 0; j < edges.size(); ++j) {
                if(edges[j].contains(points[1])) {
                    index = j;
                    break;
                }
            }

            if(index < vertices.size()) {

                for(std::size_t j = index; j < vertices.size(); ++j) {
                    if(line > vertices[j])
                        new_vertices.emplace_back(vertices[j]);
                }

                for(std::size_t j = 0; j < index; ++j) {
                    if(line > vertices[j])
                        new_vertices.emplace_back(vertices[j]);
                }

            } else {
                
                for(const auto &vertex: vertices) {
                    if(line > vertex)
                        new_vertices.emplace_back(vertex);
                }

            }

        } else { // Point above the Line.

            if((std::abs(line[0]) <= GEOMETRY_TOLERANCE) && (points[0][0] > points[1][0])) {
                std::swap(points[0], points[1]);
            } else if((line.angular() > 0.0) && (points[0][1] > points[1][1])) {
                std::swap(points[0], points[1]);
            } else if((line.angular() < 0.0) && (points[0][1] < points[1][1])) {
                std::swap(points[0], points[1]);
            }

            new_vertices.emplace_back(points[0]);
            new_vertices.emplace_back(points[1]);

            for(std::size_t j = 0; j < edges.size(); ++j) {
                if(edges[j].contains(points[1])) {
                    index = j;
                    break;
                }
            }

            if(index < vertices.size()) {

                for(std::size_t j = index; j < vertices.size(); ++j) {
                    if(line <= vertices[j])
                        new_vertices.emplace_back(vertices[j]);
                }

                for(std::size_t j = 0; j < index; ++j) {
                    if(line <= vertices[j])
                        new_vertices.emplace_back(vertices[j]);
                }

            } else {
                
                for(const auto &vertex: vertices) {
                    if(line <= vertex)
                        new_vertices.emplace_back(vertex);
                }

            }

        }

        return Polygon{new_vertices};
    }

}