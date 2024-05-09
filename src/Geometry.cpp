/**
 * @file Geometry.cpp
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

        Real px = p[0], qx = q[0];
        Real py = p[1], qy = q[1];
        Real mx = (px + qx) / 2, my = (py + qy) / 2;

        // Evaluation by cases.

        // Px = Qx.
        if(std::abs(px - qx) <= GEOMETRY_TOLERANCE)
            return Line{0.0, 1.0, my};

        // Py = Qy.
        if(std::abs(py - qy) <= GEOMETRY_TOLERANCE)
            return Line{1.0, 0.0, mx};

        // Default.
        Real a = qx - px;
        Real b = qy - py;
        Real c = a * mx + b * my;
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
        Real a = - s[0] / s[1];
        Real c = s[2] / s[1];

        // r: y = bx + d.
        Real b = -r[0] / r[1];
        Real d = r[2] / r[1];

        Real x = (d - c) / (a - b);
        Real y = a * x + c;

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

        for(const auto &point: intersections(line, segment.line()))
            if(segment.contains(point))
                points.emplace_back(point);

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

        for(const auto &segment: polygon.edges())
            for(const auto &point: intersections(line, segment))
                points.emplace_back(point);

        return points;
    }

    /**
     * @brief Polygon collapse over a vertex.
     * 
     * @param polygon 
     * @param point 
     * @return Polygon 
     */
    Polygon collapse(const Polygon &polygon, const Point &point) {
        #ifndef NDEBUG // Integrity check.
        assert(polygon.points.size() > 3);
        #endif

        std::vector<Point> vertices;

        for(const auto &vertex: polygon.points) {
            if(point == vertex)
                continue;

            vertices.emplace_back(vertex);
        }

        return Polygon{vertices};
    }

    /**
     * @brief Polygon collape over an edge.
     * 
     * @param polygon 
     * @param segment 
     * @return Polygon 
     */
    Polygon collapse(const Polygon &polygon, const Segment &segment) {
        #ifndef NDEBUG // Integrity check.
        std::size_t counter = 0;

        for(const auto &vertex: polygon.points)
            if((vertex == segment[0]) || (vertex == segment[1]))
                ++counter;

        assert(polygon.points.size() > 3);
        assert(counter == 2);
        #endif

        std::vector<Point> vertices;
        Point mid = (segment[0] + segment[1]) * 0.5;

        for(const auto &vertex: polygon.points) {
            if(vertex == segment[0])
                continue;

            if(vertex == segment[1]) {
                vertices.emplace_back(mid);
                continue;
            }

            vertices.emplace_back(vertex);
        }

        return Polygon{vertices};
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
        std::vector<Point> vertices = polygon.points;
        std::vector<Segment> edges = polygon.edges();
        std::size_t index = 0;

        if(line > point) { // Point under the Line.

            if((std::abs(line[0]) <= GEOMETRY_TOLERANCE) && (points[0][0] < points[1][0]))
                std::swap(points[0], points[1]);
            else if((line.angular() > 0.0) && (points[0][1] < points[1][1]))
                std::swap(points[0], points[1]);
            else if((line.angular() < 0.0) && (points[0][1] > points[1][1]))
                std::swap(points[0], points[1]);

            new_vertices.emplace_back(points[0]);
            new_vertices.emplace_back(points[1]);

            for(std::size_t j = 0; j < edges.size(); ++j)
                if(edges[j].contains(points[1])) {
                    index = j;
                    break;
                }

            if(index < vertices.size()) {
                for(std::size_t j = index; j < vertices.size(); ++j)
                    if(line > vertices[j])
                        new_vertices.emplace_back(vertices[j]);

                for(std::size_t j = 0; j < index; ++j)
                    if(line > vertices[j])
                        new_vertices.emplace_back(vertices[j]);

            } else
                for(const auto &vertex: vertices)
                    if(line > vertex)
                        new_vertices.emplace_back(vertex);

        } else { // Point above the Line.

            if((std::abs(line[0]) <= GEOMETRY_TOLERANCE) && (points[0][0] > points[1][0]))
                std::swap(points[0], points[1]);
            else if((line.angular() > 0.0) && (points[0][1] > points[1][1]))
                std::swap(points[0], points[1]);
            else if((line.angular() < 0.0) && (points[0][1] < points[1][1]))
                std::swap(points[0], points[1]);

            new_vertices.emplace_back(points[0]);
            new_vertices.emplace_back(points[1]);

            for(std::size_t j = 0; j < edges.size(); ++j)
                if(edges[j].contains(points[1])) {
                    index = j;
                    break;
                }

            if(index < vertices.size()) {
                for(std::size_t j = index; j < vertices.size(); ++j)
                    if(line <= vertices[j])
                        new_vertices.emplace_back(vertices[j]);

                for(std::size_t j = 0; j < index; ++j)
                    if(line <= vertices[j])
                        new_vertices.emplace_back(vertices[j]);

            } else
                for(const auto &vertex: vertices)
                    if(line <= vertex)
                        new_vertices.emplace_back(vertex);

        }

        return Polygon{new_vertices};
    }

}

namespace std {

    /**
     * @brief Return the norm of a Point.
     * 
     * @param point 
     * @return Real 
     */
    pacs::Real abs(const pacs::Point &point) {
        return std::sqrt(point[0] * point[0] + point[1] * point[1]);
    }

    /**
     * @brief Return the length of a Segment.
     * 
     * @param segment 
     * @return Real 
     */
    pacs::Real abs(const pacs::Segment &segment) {
        return std::abs(segment[1] - segment[0]);
    }

}