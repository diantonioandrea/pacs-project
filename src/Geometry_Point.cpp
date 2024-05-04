/**
 * @file Geometry_Point.cpp
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

}