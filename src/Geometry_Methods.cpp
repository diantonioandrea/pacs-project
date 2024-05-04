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
        assert(!(p - q).is_origin());
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
        double m = (px - qx) / (py - qy);
        return Line{m, 1.0, m * mx + my};
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
        double x = (s[2] / s[1] - r[2] / r[1]) / (s[0] / s[1] - r[0] / r[1]);
        double y = s[2] / s[1] - s[0] / s[1] * x;

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

        #ifndef NDEBUG // Integrity check.
        assert(polygon.contains(point));
        assert((points.size() == 0) || (points.size() == 2));
        #endif

        if(points.size() == 0)
            return polygon;
        
        // Convexity assumption.
        std::vector<Point> lower{points}, upper{points};

        for(const auto &point: polygon.vertices()) {
            if(point < line) {
                lower.emplace_back(point);
            } else {
                upper.emplace_back(point);
            }
        }

        // Convexity sorting.
        Point lower_centroid = Polygon{lower}.centroid();
        Point upper_centroid = Polygon{upper}.centroid();

        for(std::size_t j = 0; j < lower.size(); ++j) {
            for(std::size_t k = 0; k < lower.size() - j - 1; ++k) {
                if(Segment{lower_centroid, lower[k]}.orientation() > Segment{lower_centroid, lower[k + 1]}.orientation()) {
                    Point temp = lower[k];
                    lower[k] = lower[k + 1]; 
                    lower[k + 1] = temp;
                }
            }
        }

        for(std::size_t j = 0; j < upper.size(); ++j) {
            for(std::size_t k = 0; k < upper.size() - j - 1; ++k) {
                if(Segment{upper_centroid, upper[k]}.orientation() > Segment{upper_centroid, upper[k + 1]}.orientation()) {
                    Point temp = upper[k];
                    upper[k] = upper[k + 1]; 
                    upper[k + 1] = temp;
                }
            }
        }

        // Interior check.
        if(Polygon{lower}.contains(point))
            return Polygon{lower};

        return Polygon{upper};
    }

}