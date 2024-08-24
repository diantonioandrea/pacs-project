/**
 * @file Geometry.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <cassert>
#include <cmath>
#include <algorithm>

namespace pacs {
    
    /**
     * @brief Returns the distance between two Points.
     * 
     * @param p 
     * @param q 
     * @return Real 
     */
    Real distance(const Point &p, const Point &q) {
        if(p == q)
            return 0.0;
            
        return std::abs(p - q);
    }

    /**
     * @brief Returns the distance between a Point and a Line.
     * 
     * @param point 
     * @param line 
     * @return Real 
     */
    Real distance(const Point &point, const Line &line) {
        if(line.contains(point))
            return 0.0;

        // Evaluation by cases.
        
        // Horizonal line.
        if(std::abs(line[0]) < GEOMETRY_TOLERANCE)
            return std::abs(point[1] - line[2] / line[1]);

        // Vertical line.
        if(std::abs(line[1]) < GEOMETRY_TOLERANCE)
            return std::abs(point[0] - line[2] / line[0]);

        // General case.
        Point intersection = intersections(line, normal(line, point))[0];

        return distance(point, intersection);
    }

    Real distance(const Point &point, const Segment &segment) {
        if(segment.contains(point))
            return 0.0;

        Real d0 = distance(point, segment[0]), d1 = distance(point, segment[1]);

        // General case.
        Line line = segment.line();
        Point intersection = intersections(line, normal(segment, point))[0];

        if(segment.contains(intersection))
            return distance(point, intersection);
        
        return (d0 < d1) ? d0 : d1;
    }
    
    /**
     * @brief Finds the bisector given two Points.
     * 
     * @param p 
     * @param q 
     * @return Line 
     */
    Line bisector(const Point &p, const Point &q) {
        #ifndef NDEBUG // Checks.
        assert(p != q);
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

        // No (or infinite) intersections due to parallelism.
        if((std::abs(s[0] - r[0]) <= GEOMETRY_TOLERANCE) && (std::abs(s[1] - r[1]) <= GEOMETRY_TOLERANCE))
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

        // S is horizontal.
        if(std::abs(s[0]) <= GEOMETRY_TOLERANCE) {
            points.emplace_back(r.x(s[2] / s[1]), s[2] / s[1]);
            return points;
        }

        // R is horizontal.
        if(std::abs(r[0]) <= GEOMETRY_TOLERANCE) {
            points.emplace_back(s.x(r[2] / r[1]), r[2] / r[1]);
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
        std::vector<Point> candidates;
        std::vector<Point> points;

        for(const auto &segment: polygon.edges())
            for(const auto &point: intersections(line, segment))
                candidates.emplace_back(point);

        for(std::size_t j = 0; j < candidates.size(); ++j) {
            bool keep = true;

            for(std::size_t k = 0; k < points.size(); ++k) {
                if(candidates[j] == points[k]) {
                    keep = false;
                    break;
                }
            }

            if(keep)
                points.emplace_back(candidates[j]);
        }

        return points;
    }

    /**
     * @brief Returns the normal with respect to a Line and passing through a point.
     * 
     * @param line 
     * @param point 
     * @return Line 
     */
    Line normal(const Line &line, const Point &point) {
        // Evaluation by cases.

        // Horizonal line.
        if(std::abs(line[0]) < GEOMETRY_TOLERANCE)
            return Line{1.0, 0.0, point[0]};

        // Vertical line.
        if(std::abs(line[1]) < GEOMETRY_TOLERANCE)
            return Line{0.0, 1.0, point[1]};

        // General case.
        Real m = line[1] / line[0];
        return Line{-m, 1.0, -m * point[0] + point[1]};
    }

    /**
     * @brief Returns the normal with respect to a Segment and passing through a point.
     * 
     * @param segment 
     * @param point 
     * @return Line 
     */
    Line normal(const Segment &segment, const Point &point) {
        return normal(segment.line(), point);
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
     * @brief Returns the reflection(s) of a Point with respect to a Polygon.
     * 
     * @param polygon 
     * @param point 
     * @return Point 
     */
    std::vector<Point> reflections(const Polygon &polygon, const Point &point) {
        #ifndef NDEBUG // Integrity check.
        assert(polygon.contains(point));
        #endif
        
        std::vector<Point> points;

        // Reflections with respect to the edges.
        std::vector<Segment> edges = polygon.edges();

        for(std::size_t j = 0; j < polygon.points.size(); ++j) {
            std::vector<Point> normal_intersections = intersections(normal(edges[j], point), edges[j]);

            if(normal_intersections.size() == 0)
                continue;

            Point reflection = normal_intersections[0] * 2 - point;

            if(!(polygon.contains(reflection)))
                points.emplace_back(reflection);
        }

        // Reflection with respect to the vertices.
        for(std::size_t j = 0; j < polygon.points.size(); ++j) {
            Point reflection = polygon.points[j] * 2 - point;

            if(!(polygon.contains(reflection)))
                points.emplace_back(reflection);
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
        std::vector<Point> candidates = intersections(line, polygon);
        std::vector<Point> new_vertices;

        #ifndef NDEBUG // Integrity check.
        assert(polygon.contains(point));
        #endif

        if(candidates.size() <= 1)
            return polygon;

        // Building.
        std::vector<Point> vertices = polygon.points;
        std::vector<Segment> edges = polygon.edges();
        std::array<std::size_t, 2> indices{0, 1};

        // Points.
        std::array<Point, 2> points{candidates[0], candidates[1]};

        if(candidates.size() > 2) // Handles candidates.size() >= 3.
            for(std::size_t j = 0; j < candidates.size(); ++j)
                for(std::size_t k = 0; k < candidates.size(); ++k) {
                    if(j == k)
                        continue;

                    if(std::abs(candidates[j] - candidates[k]) > std::abs(points[0] - points[1]))
                        points = {candidates[j], candidates[k]};
                }

        // Indices.
        for(std::size_t j = 0; j < edges.size(); ++j) {
            Segment edge = edges[j];

            if(edge.contains(points[0]) && edge.contains(points[1])) // Shouldn't happen.
                return polygon;

            if(edge.contains(points[0]) && (edge[1] != points[0]))
                indices[0] = j;
            
            if(edge.contains(points[1]) && (edge[1] != points[1]))
                indices[1] = j;
        }

        if(indices[0] > indices[1]) {
            std::swap(indices[0], indices[1]);
            std::swap(points[0], points[1]);
        }

        // New polygons.
        std::vector<Point> a_points, b_points;

        for(std::size_t j = 0; j < vertices.size(); ++j) {
            if((j <= indices[0]) || (j > indices[1]))
                a_points.emplace_back(vertices[j]);

            if(j == indices[0]) {
                if(points[0] != *--a_points.end()) // Avoids duplicates.
                    a_points.emplace_back(points[0]);

                b_points.emplace_back(points[0]);
            }

            if((j > indices[0]) && (j <= indices[1]))
                if(vertices[j] != *--b_points.end()) // Avoids duplicates.
                    b_points.emplace_back(vertices[j]);

            if(j == indices[1]) {
                a_points.emplace_back(points[1]);
                b_points.emplace_back(points[1]);
            }
        }

        Polygon a_polygon{a_points}, b_polygon{b_points};

        if(a_polygon.contains(point))
            return a_polygon;

        #ifndef NDEBUG // Integrity check.
        assert(b_polygon.contains(point));
        #endif

        return b_polygon;
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