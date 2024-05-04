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

    // READ.

    /**
     * @brief Returns the j-th parameter.
     * 
     * @param j 
     * @return double 
     */
    double Line::operator [](const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j <= 2);
        #endif

        return (j == 0) ? this->a : ((j == 1) ? this->b : this->c);
    }

    /**
     * @brief Returns the first coordinate of a Line's point given the second one.
     * 
     * @param y 
     * @return double 
     */
    double Line::x(const double &y) const {
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(this->a) > GEOMETRY_TOLERANCE);
        #endif

        return (this->c - this->b * y) / this->a;
    }

    /**
     * @brief Returns the second coordinate of a Line's point given the first one.
     * 
     * @param x 
     * @return double 
     */
    double Line::y(const double &x) const {
        #ifndef NDEBUG // Integrity check.
        assert(std::abs(this->b) > GEOMETRY_TOLERANCE);
        #endif

        return (this->c - this->a * x) / this->b;
    }

    // METHODS.

    /**
     * @brief Checks whether a Point is contained inside the Line.
     * 
     * @param point 
     * @return true 
     * @return false 
     */
    bool Line::contains(const Point &point) const {
        return std::abs(this->a * point[0] + this->b * point[1] - this->c) <= GEOMETRY_TOLERANCE;
    }

    /**
     * @brief Checks whether two Lines are parallel.
     * 
     * @param line 
     * @return true 
     * @return false 
     */
    bool Line::is_parallel(const Line &line) const {
        return (std::abs(this->a - line.a) <= GEOMETRY_TOLERANCE) && (std::abs(this->b - line.b) <= GEOMETRY_TOLERANCE);
    }

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

    // OUTPUT.

    std::ostream &operator <<(std::ostream &ost, const Segment &segment) {
        return ost << "[" << segment.a << ", " << segment.b << "]";
    }


    // POLYGON.

    // CONSTRUCTORS.
    
    /**
     * @brief Constructs a new Polygon from a given vector of points.
     * 
     * @param points 
     */
    Polygon::Polygon(const std::vector<Point> &points): points{points} {}

    /**
     * @brief Copy constructor.
     * 
     * @param polygon 
     */
    Polygon::Polygon(const Polygon &polygon): points{polygon.points} {}

    // METHODS.
    
    /**
     * @brief Returns the vector of Segments.
     * 
     * @return std::vector<Segment> 
     */
    std::vector<Segment> Polygon::segments() const {
        std::vector<Segment> segments;

        for(std::size_t j = 0; j < this->points.size() - 1; ++j)
            segments.emplace_back(this->points[j], this->points[j + 1]);

        segments.emplace_back(*--(this->points.end()), this->points[0]);

        return segments;
    }

    /**
     * @brief Checks whether a point is contained inside the Polygon. Does not count Points over the perimeter.
     * 
     * @param point 
     * @return true 
     * @return false 
     */
    bool Polygon::contains(const Point &point) const {
        for(const auto &segment: this->segments()) {
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

        for(const auto &segment: polygon.segments()) {
            for(const auto &point: intersections(line, segment))
                points.emplace_back(point);
        }

        return points;
    }

}