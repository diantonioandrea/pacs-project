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

    std::vector<Point> Polygon::vertices() const {
        return this->points;
    }
    
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
        std::vector<Point> points;

        for(const auto &intersection: intersections(horizonal, *this)) {
            if(intersection[0] > point[0])
                points.emplace_back(intersection);
        }

        return points.size() % 2 == 1;
    }

    /**
     * @brief Returns the polygon's centroid.
     * 
     * @return Point 
     */
    Point Polygon::centroid() const {
        Point sum{0.0, 0.0};

        for(const auto &point: this->points)
            sum += point;

        return sum *= (1.0 / this->points.size());
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
        ost << "{";

        for(std::size_t j = 0; j < polygon.points.size(); ++j) {
            ost << polygon.points[j];

            if(j < polygon.points.size() - 1)
                ost << ", ";
        }

        return ost;
    }

}