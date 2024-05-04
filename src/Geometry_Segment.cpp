/**
 * @file Geometry_Segment.cpp
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
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Segment from a given a and b.
     * 
     * @param a 
     * @param b 
     */
    Segment::Segment(const Point &a, const Point &b): a{a}, b{b} {
        #ifndef NDEBUG // Integrity check.
        assert(!(a - b).is_origin());
        #endif
    }

    /**
     * @brief Constructs a new Segment from a given array of points.
     * 
     * @param ab 
     */
    Segment::Segment(const std::array<Point, 2> &ab): a{ab[0]}, b{ab[1]} {
        #ifndef NDEBUG // Integrity check.
        assert(!(ab[0] - ab[0]).is_origin());
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param segment 
     */
    Segment::Segment(const Segment &segment): a{segment.a}, b{segment.b} {}

    // METHODS.

    /**
     * @brief Returs the Line passing through the Segment's extremes.
     * 
     * @return Line 
     */
    Line Segment::line() const {
        double ax = this->a[0], bx = this->b[0];
        double ay = this->a[1], by = this->b[1];

        // Evaluation by cases.

        // Ax = Bx.
        if(std::abs(ax - bx) <= GEOMETRY_TOLERANCE)
            return Line{1.0, 0.0, ax};

        // Ay = By.
        if(std::abs(ay - by) <= GEOMETRY_TOLERANCE)
            return Line{0.0, 1.0, ay};

        // Default.
        return Line{(ay - by) / (bx - ax), 1.0, (ay - by) / (bx - ax) * ax + ay};
    }

    /**
     * @brief Checks whether a Point is contained inside the Segment.
     * 
     * @param point 
     * @return true 
     * @return false 
     */
    bool Segment::contains(const Point &point) const {
        if((this->a - point).is_origin() || (this->b - point).is_origin())
            return true;

        if(!(this->line().contains(point)))
            return false;

        bool x = ((this->a[0] <= point[0]) && (point[0] <= this->b[0])) || ((this->b[0] <= point[0]) && (point[0] <= this->a[0]));
        bool y = ((this->a[1] <= point[1]) && (point[1] <= this->b[1])) || ((this->b[1] <= point[1]) && (point[1] <= this->a[1]));

        return x && y;
    }

    /**
     * @brief Returns the Segment's orientation.
     * 
     * @return double 
     */
    double Segment::orientation() const {
        Point vector = this->b - this->a;

        // Vertical vector.
        if(std::abs(vector[0]) <= GEOMETRY_TOLERANCE)
            return (vector[1] >= 0) ? M_PI_2 : 1.5 * M_PI;

        double angle = std::atan(vector[1] / vector[0]);

        if(vector[1] >= 0) {
            if(vector[0] < 0)
                angle += M_PI;
        } else {
            if(vector[0] >= 0)
                angle += 2.0 * M_PI;
            else
                angle += M_PI;
        }

        return angle;
    }

    // OUTPUT.

    std::ostream &operator <<(std::ostream &ost, const Segment &segment) {
        return ost << "[" << segment.a << ", " << segment.b << "]";
    }

}