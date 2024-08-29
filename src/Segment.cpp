/**
 * @file Segment.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-05-04
 *
 * @copyright Copyright (c) 2024
 *
 */
// clang-format off
/*
Handling geometric predicates is tricky and sometimes one wants to use the maximal precision to avoid issues.
To this aim, some libraries may be of help such as:
- CGAL (https://www.cgal.org/)
- Shewchuk's predicates (http://www.cs.cmu.edu/~quake/robust.html)
- Predicates for computational geometry (http://www.cs.cmu.edu/~quake/robust.html)
*/
#include <PacsHPDG.hpp>

#include <cassert>
#include <cmath>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Segment from a given a and b.
     * 
     * @param a Point.
     * @param b Point.
     */
    Segment::Segment(const Point &a, const Point &b): a{a}, b{b} {
        #ifndef NDEBUG // Integrity check.
        assert(a != b);
        #endif
    }

    /**
     * @brief Constructs a new Segment from a given array of points.
     * 
     * @param ab Points.
     */
    Segment::Segment(const std::array<Point, 2> &ab): a{ab[0]}, b{ab[1]} {
        #ifndef NDEBUG // Integrity check.
        assert(ab[0] != ab[1]);
        #endif
    }

    /**
     * @brief Copy constructor.
     * 
     * @param segment Segment.
     */
    Segment::Segment(const Segment &segment): a{segment.a}, b{segment.b} {}

    // READ.

    /**
     * @brief Returns the j-th extreme.
     * 
     * @param j Index.
     * @return Point 
     */
    Point Segment::operator [](const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        return (j == 0) ? this->a : this->b;
    }

    // COMPARISONS

    /**
     * @brief Segment == Segment.
     * 
     * @param segment Segment.
     * @return true 
     * @return false 
     */
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
        Real ax = this->a[0], bx = this->b[0];
        Real ay = this->a[1], by = this->b[1];

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
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Segment::contains(const Point &point) const {
        if((this->a == point) || (this->b == point))
            return true;

        if(!(this->line().contains(point)))
            return false;

        // Vertical segment.
        if(std::abs(this->a[0] - this->b[0]) <= GEOMETRY_TOLERANCE) {
            if(std::abs(this->a[0] - point[0]) > GEOMETRY_TOLERANCE)
                return false;

            if(this->a[1] <= this->b[1])
                return (this->a[1] <= point[1]) && (point[1] <= this->b[1]);
            
            return (this->b[1] <= point[1]) && (point[1] <= this->a[1]);
        }

        // Horizontal segment.
        if(std::abs(this->a[1] - this->b[1]) <= GEOMETRY_TOLERANCE) {
            if(std::abs(this->a[1] - point[1]) > GEOMETRY_TOLERANCE)
                return false;

            if(this->a[0] <= this->b[0])
                return (this->a[0] <= point[0]) && (point[0] <= this->b[0]);
            
            return (this->b[0] <= point[0]) && (point[0] <= this->a[0]);
        }

        // General cases.
        if(this->a[0] <= this->b[0]) {
            if(this->a[1] <= this->b[1])
                return (this->a[0] <= point[0]) && (point[0] <= this->b[0]) && (this->a[1] <= point[1]) && (point[1] <= this->b[1]);

            return (this->a[0] <= point[0]) && (point[0] <= this->b[0]) && (this->b[1] <= point[1]) && (point[1] <= this->a[1]);
        }

        if(this->a[1] <= this->b[1])
            return (this->b[0] <= point[0]) && (point[0] <= this->a[0]) && (this->a[1] <= point[1]) && (point[1] <= this->b[1]);

        return (this->b[0] <= point[0]) && (point[0] <= this->a[0]) && (this->b[1] <= point[1]) && (point[1] <= this->a[1]);
    }

    /**
     * @brief Checks whether a Segment is contained inside the Segment.
     * 
     * @param segment Segment.
     * @return true 
     * @return false 
     */
    bool Segment::contains(const Segment &segment) const {
        return this->contains(segment.a) && this->contains(segment.b);
    }

    // OUTPUT.

    /**
     * @brief Segment output.
     * 
     * @param ost 
     * @param segment 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Segment &segment) {
        return ost << "[" << segment.a << ", " << segment.b << "]";
    }

}