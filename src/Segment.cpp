/**
 * @file Segment.cpp
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
        assert(!(a - b).is_zero());
        #endif
    }

    /**
     * @brief Constructs a new Segment from a given array of points.
     * 
     * @param ab 
     */
    Segment::Segment(const std::array<Point, 2> &ab): a{ab[0]}, b{ab[1]} {
        #ifndef NDEBUG // Integrity check.
        assert(!(ab[0] - ab[0]).is_zero());
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param segment 
     */
    Segment::Segment(const Segment &segment): a{segment.a}, b{segment.b} {}

    // READ.

    /**
     * @brief Returns the j-th extreme.
     * 
     * @param j 
     * @return Point 
     */
    Point Segment::operator [](const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        return (j == 0) ? this->a : this->b;
    }

    // COMPARISONS


    bool Segment::operator ==(const Segment &segment) const {
        return ((this->a == segment.a) && (this->b == segment.b)) || ((this->a == segment.b) && (this->b == segment.a));
    }

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
        if((this->a - point).is_zero() || (this->b - point).is_zero())
            return true;

        if(!(this->line().contains(point)))
            return false;

        double min_x = this->a[0] <= this->b[0] ? this->a[0] : this->b[0];
        double max_x = this->a[0] >= this->b[0] ? this->a[0] : this->b[0];

        double min_y = this->a[1] <= this->b[1] ? this->a[1] : this->b[1];
        double max_y = this->a[1] >= this->b[1] ? this->a[1] : this->b[1];

        bool x = (min_x <= point[0]) && (point[0] <= max_x);
        bool y = (min_y <= point[1]) && (point[1] <= max_y);

        if(std::abs(this->a[0] - this->b[0]) <= GEOMETRY_TOLERANCE)
            x = std::abs(this->a[0] - point[0]) <= GEOMETRY_TOLERANCE;

        if(std::abs(this->a[1] - this->b[1]) <= GEOMETRY_TOLERANCE)
            y = std::abs(this->a[1] - point[1]) <= GEOMETRY_TOLERANCE;

        return x && y;
    }

    bool Segment::contains(const Segment &segment) const {
        return this->contains(segment.a) && this->contains(segment.b);
    }

    // OUTPUT.

    std::ostream &operator <<(std::ostream &ost, const Segment &segment) {
        return ost << "[" << segment.a << ", " << segment.b << "]";
    }

}