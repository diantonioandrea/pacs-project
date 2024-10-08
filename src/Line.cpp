/**
 * @file Line.cpp
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

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Line from given a, b and c.
     * 
     * @param a Parameter.
     * @param b Parameter.
     * @param c Parameter.
     */
    Line::Line(const Real &a, const Real &b, const Real &c): a{a}, b{b}, c{c} {}

    /**
     * @brief Constructs a new Line from an array of parameters.
     * 
     * @param abc Parameters.
     */
    Line::Line(const std::array<Real, 3> &abc): a{abc[0]}, b{abc[1]}, c{abc[2]} {}

    /**
     * @brief Copy constructor.
     * 
     * @param line Line.
     */
    Line::Line(const Line &line): a{line.a}, b{line.b}, c{line.c} {}

    // READ.

    /**
     * @brief Returns the j-th parameter.
     * 
     * @param j Index.
     * @return Real 
     */
    Real Line::operator [](const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        #endif

        return (j == 0) ? this->a : ((j == 1) ? this->b : this->c);
    }

    /**
     * @brief Returns the first coordinate of a Line's point given the second one.
     * 
     * @param y Coordinate.
     * @return Real 
     */
    Real Line::x(const Real &y) const {
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(this->a) > GEOMETRY_TOLERANCE);
        #endif

        return (this->c - this->b * y) / this->a;
    }

    /**
     * @brief Returns the second coordinate of a Line's point given the first one.
     * 
     * @param x Coordinate.
     * @return Real 
     */
    Real Line::y(const Real &x) const {
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(this->b) > GEOMETRY_TOLERANCE);
        #endif

        return (this->c - this->a * x) / this->b;
    }

    // COMPARISONS.

    /**
     * @brief Line == Line.
     * 
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool Line::operator ==(const Line &line) const {
        bool a = std::abs(this->a - line.a) <= GEOMETRY_TOLERANCE;
        bool b = std::abs(this->b - line.b) <= GEOMETRY_TOLERANCE;
        bool c = std::abs(this->c - line.c) <= GEOMETRY_TOLERANCE;

        return a && b && c;
    }

    /**
     * @brief Line < Point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Line::operator <(const Point &point) const {

        // Line is vertical.
        if(std::abs(this->b) <= GEOMETRY_TOLERANCE)
            return this->c / this->a < point[0];

        // Line is horizontal.
        if(std::abs(this->a) <= GEOMETRY_TOLERANCE)
            return this->c / this->b < point[1];

        // Vertical Line intersection.
        Line vertical{1.0, 0.0, point[0]};
        Point intersection = intersections(*this, vertical)[0];

        return intersection[1] < point[1];
    }

    /**
     * @brief Line <= Point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Line::operator <=(const Point &point) const {
        return (*this < point) || (this->contains(point));
    }

    /**
     * @brief Line > Point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Line::operator >(const Point &point) const {
        return !(*this <= point);
    }

    /**
     * @brief Line >= Point.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Line::operator >=(const Point &point) const {
        return !(*this < point);
    }

    /**
     * @brief Point < Line.
     * 
     * @param point Point.
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool operator <(const Point &point, const Line &line) {
        return line > point;
    }

    /**
     * @brief Point <= Line.
     * 
     * @param point Point.
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool operator <=(const Point &point, const Line &line) {
        return line >= point;
    }

    /**
     * @brief Point > Line.
     * 
     * @param point Point.
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool operator >(const Point &point, const Line &line) {
        return line < point;
    }

    /**
     * @brief Point >= Line.
     * 
     * @param point Point.
     * @param line Line.
     * @return true 
     * @return false 
     */
    bool operator >=(const Point &point, const Line &line) {
        return line <= point;
    }

    // METHODS.

    /**
     * @brief Checks whether a Point is contained inside the Line.
     * 
     * @param point Point.
     * @return true 
     * @return false 
     */
    bool Line::contains(const Point &point) const {
        // Horizontal line.
        if(std::abs(this->a) <= GEOMETRY_TOLERANCE)
            return std::abs(point[1] - this->c / this->b) <= GEOMETRY_TOLERANCE;

        // Vertical line.
        if(std::abs(this->b) <= GEOMETRY_TOLERANCE)
            return std::abs(point[0] - this->c / this->a) <= GEOMETRY_TOLERANCE;

        // General case.
        return std::abs(this->a * point[0] + this->b * point[1] - this->c) <= GEOMETRY_TOLERANCE;
    }

    // OUTPUT.

    /**
     * @brief Line output.
     * 
     * @param ost 
     * @param line Line.
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Line &line) {
        return ost << "(" << line.a << ") x + " << "(" << line.b << ") y = " << "(" << line.c << ")";
    }
    
}