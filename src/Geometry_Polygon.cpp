/**
 * @file Geometry_Polygon.cpp
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

// Random.
#include <cstdlib>
#include <ctime>

namespace pacs {

    // CONSTRUCTORS.
    
    /**
     * @brief Constructs a new Polygon from a given vector of points.
     * 
     * @param points 
     */
    Polygon::Polygon(const std::vector<Point> &points): points{points} {
        #ifndef NDEBUG // Integrity check.
        assert(points.size() > 2);
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param polygon 
     */
    Polygon::Polygon(const Polygon &polygon): points{polygon.points} {}

    /**
     * @brief Copy operator.
     * 
     * @param polygon 
     * @return Polygon& 
     */
    Polygon &Polygon::operator =(const Polygon &polygon) {
        this->points = polygon.points;

        return *this;
    }

    // METHODS.
    
    /**
     * @brief Returns the vector of Segments.
     * 
     * @return std::vector<Segment> 
     */
    std::vector<Segment> Polygon::edges() const {
        std::vector<Segment> edges;

        for(std::size_t j = 0; j < this->points.size() - 1; ++j)
            edges.emplace_back(this->points[j], this->points[j + 1]);

        edges.emplace_back(*--(this->points.end()), this->points[0]);

        return edges;
    }
    
    /**
     * @brief Checks whether a point is contained inside the Polygon. Does not count Points over the perimeter.
     * 
     * @param point 
     * @return true 
     * @return false 
     */
    bool Polygon::contains(const Point &point) const {
        for(const auto &segment: this->edges()) {
            if(segment.contains(point))
                return false;
        }

        // Horizontal line.
        Line horizonal{0.0, 1.0, point[1]};
        std::size_t points = 0;

        for(const auto &intersection: intersections(horizonal, *this)) {
            if(intersection[0] > point[0])
                points += 1;
        }

        return points % 2;
    }

    /**
     * @brief Checks whether a segment is contained inside the Polygon's perimeter.
     * 
     * @return true 
     * @return false 
     */
    bool Polygon::contains(const Segment &segment) const {
        for(const auto &edge: this->edges()) {
            if(edge == segment)
                return true;
        }

        return false;
    }

    /**
     * @brief Returns the Polygon's area.
     * 
     * @return double 
     */
    double Polygon::area() const {
        double area = 0.0;

        for(const auto &edge: this->edges())
            area += (edge[0][1] + edge[1][1]) * (edge[0][0] - edge[1][0]);

        return 0.5 * area;
    }

    /**
     * @brief Returns the Polygon's centroid.
     * 
     * @return Point 
     */
    Point Polygon::centroid() const {
        Point centroid{0.0, 0.0};

        for(const auto &edge: this->edges()) {
            centroid[0] += (edge[0][0] + edge[1][0]) * (edge[0][0] * edge[1][1] - edge[1][0] * edge[0][1]);
            centroid[1] += (edge[0][1] + edge[1][1]) * (edge[0][0] * edge[1][1] - edge[1][0] * edge[0][1]);
        }

        return centroid * (1 / (6 * this->area()));
    }

    /**
     * @brief Returns a random point inside the Polygon.
     * 
     * @return Point 
     */
    Point Polygon::random() const {
        // Boundaries.
        double x_min = this->points[0][0], x_max = this->points[0][0];
        double y_min = this->points[0][1], y_max = this->points[0][1];
        double x, y;

        for(std::size_t j = 1; j < this->points.size(); ++j) {
            x_min = (points[j][0] < x_min) ? points[j][0] : x_min;
            x_max = (points[j][0] > x_max) ? points[j][0] : x_max;
            y_min = (points[j][1] < y_min) ? points[j][1] : y_min;
            y_max = (points[j][1] > y_max) ? points[j][1] : y_max;
        }

        // Generation.
        do {
            x = x_min + (x_max - x_min) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
            y = y_min + (y_max - y_min) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        } while(!this->contains(Point{x, y}));

        return Point{x, y};
    }

    // OUTPUT.

    /**
     * @brief Polygon output.
     * 
     * @param ost 
     * @param polygon 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Polygon &polygon) {
        ost << "[";

        for(std::size_t j = 0; j < polygon.points.size(); ++j) {
            ost << polygon.points[j];

            if(j < polygon.points.size() - 1)
                ost << ", ";
        }

        return ost << "]";
    }

}