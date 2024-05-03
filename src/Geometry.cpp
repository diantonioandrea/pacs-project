/**
 * @file Geometry.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
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

    // POINT.

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Point from given x and y.
     * 
     * @param x 
     * @param y 
     */
    Point::Point(const double &x, const double &y): x{x}, y{y} {}

    /**
     * @brief Constructs a new Point from an array of coordinates.
     * 
     * @param xy 
     */
    Point::Point(const std::array<double, 2> &xy): x{xy[0]}, y{xy[1]} {}

    /**
     * @brief Copy constructor.
     * 
     * @param point 
     */
    Point::Point(const Point &point): x{point.x}, y{point.y} {}

    /**
     * @brief Copy operator.
     * 
     * @param point 
     * @return Point& 
     */
    Point &Point::operator =(const Point &point) {
        this->x = point.x;
        this->y = point.y;

        return *this;
    }

    // READ AND WRITE.

    /**
     * @brief Const subscript operator.
     * 
     * @param j 
     * @return double 
     */
    double Point::operator [](const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        return (j == 0) ? this->x : this->y;
    }

    /**
     * @brief Subscript operator.
     * 
     * @param j 
     * @return double 
     */
    double &Point::operator [](const std::size_t &j) {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 1);
        #endif

        return (j == 0) ? this->x : this->y;
    }

    // METHODS.

    double Point::norm() const {
        return std::sqrt(this->x * this->x + this->y * this->y);
    }

    bool Point::is_origin() const {
        return this->norm() <= GEOMETRY_TOLERANCE;
    }

    // OPERATORS.

    /**
     * @brief Point scalar product.
     * 
     * @param scalar 
     * @return Point 
     */
    Point Point::operator *(const double &scalar) const {
        return Point{this->x * scalar, this->y * scalar};
    }

    /**
     * @brief Friend Point scalar product.
     * 
     * @param scalar 
     * @param point 
     * @return Point 
     */
    Point operator *(const double &scalar, const Point &point) {
        return Point{point.x * scalar, point.y * scalar};
    }

    /**
     * @brief Point scalar product and assignation.
     * 
     * @param scalar 
     * @return Point& 
     */
    Point &Point::operator *=(const double &scalar) {
        this->x *= scalar;
        this->y *= scalar;

        return *this;
    }

    /**
     * @brief Point sum.
     * 
     * @param point 
     * @return Point 
     */
    Point Point::operator +(const Point &point) const {
        return Point{this->x + point.x, this->y + point.y};
    }

    /**
     * @brief Point sum and assignation.
     * 
     * @param point 
     * @return Point& 
     */
    Point &Point::operator +=(const Point &point) {
        this->x += point.x;
        this->y += point.y;

        return *this;
    }

    /**
     * @brief Point difference.
     * 
     * @param point 
     * @return Point 
     */
    Point Point::operator -(const Point &point) const {
        return Point{this->x - point.x, this->y - point.y};
    }

    /**
     * @brief Point difference and assignation.
     * 
     * @param point 
     * @return Point& 
     */
    Point &Point::operator -=(const Point &point) {
        this->x -= point.x;
        this->y -= point.y;

        return *this;
    }

    // OUTPUT.

    /**
     * @brief Point output.
     * 
     * @param ost 
     * @param point 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Point &point) {
        return ost << "(" << point.x << ", " << point.y << ")";
    }


    // LINE.

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Line from given a, b and c.
     * 
     * @param a 
     * @param b 
     * @param c 
     */
    Line::Line(const double &a, const double &b, const double &c): a{a}, b{b}, c{c} {}

    /**
     * @brief Constructs a new Line from an array of parameters.
     * 
     * @param abc 
     */
    Line::Line(const std::array<double, 3> &abc): a{abc[0]}, b{abc[1]}, c{abc[2]} {}

    /**
     * @brief Copy constructor.
     * 
     * @param line 
     */
    Line::Line(const Line &line): a{line.a}, b{line.b}, c{line.c} {}

    // OUTPUT.

    /**
     * @brief Line output.
     * 
     * @param ost 
     * @param line 
     * @return std::ostream& 
     */
    std::ostream &operator <<(std::ostream &ost, const Line &line) {
        return ost << "(" << line.a << ") x + " << "(" << line.b << ") y = " << "(" << line.c << ")";
    }

    // SEGMENT.

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
     * @brief Returs the line passing through the Segment's extremes.
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

    // OUTPUT.

    std::ostream &operator <<(std::ostream &ost, const Segment &segment) {
        return ost << "[" << segment.a << ", " << segment.b << "]";
    }


    // METHODS

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

        Point mid = (p + q) * 0.5;

        double mx = mid[0], my = mid[1];
        double px = p[0], qx = q[0];
        double py = p[1], qy = q[1];

        // Evaluation by cases.

        // Px = Qx.
        if(std::abs(px - qx) <= GEOMETRY_TOLERANCE)
            return Line{0.0, 1.0, my};

        // Py = Qy.
        if(std::abs(py - qy) <= GEOMETRY_TOLERANCE)
            return Line{1.0, 0.0, mx};

        // Default.
        return Line{(px - qx) / (py - qy), 1.0, (px - qx) / (py - qy) * mx + my};
    }

}