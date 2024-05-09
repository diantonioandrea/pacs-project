/**
 * @file Geometry.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY_PACS
#define GEOMETRY_PACS

// Output.
#include <iostream>

// Containers.
#include <vector>
#include <array>

// Geometry tolerance.
#ifndef GEOMETRY_TOLERANCE
#define GEOMETRY_TOLERANCE 1E-8
#endif

namespace pacs {
    
    // Point, implemented under src/Point.cpp

    /**
     * @brief 2D Point class.
     * 
     */
    class Point {
        protected:

            // Coordinates.
            double x;
            double y;

        public:

            // CONSTRUCTORS.

            Point(const double &, const double &);
            Point(const std::array<double, 2> &);
            Point(const Point &);
            Point &operator =(const Point &);

            // READ AND WRITE.

            double operator [](const std::size_t &) const;
            double &operator [](const std::size_t &);

            // COMPARISONS.

            bool operator ==(const Point &) const;
            bool operator !=(const Point &) const;

            // METHODS.

            bool is_zero() const;

            // OPERATORS.

            Point operator *(const double &) const;
            friend Point operator *(const double &, const Point &);
            Point &operator *=(const double &);
            Point operator +(const Point &) const;
            Point &operator +=(const Point &);
            Point operator -(const Point &) const;
            Point &operator -=(const Point &);

            // OUTPUT.

            friend std::ostream &operator <<(std::ostream &, const Point &);
    };

    // Line, implemented under src/Line.cpp

    /**
     * @brief 2D Line (ax + by = c) class.
     * 
     */
    class Line {
        protected:

            // Parameters.
            const double a;
            const double b;
            const double c;

        public:

            // CONSTRUCTORS.

            Line(const double &, const double &, const double &);
            Line(const std::array<double, 3> &);
            Line(const Line &);

            // READ.

            double operator [](const std::size_t &) const;
            double x(const double &) const;
            double y(const double &) const;

            // COMPARISONS.

            bool operator ==(const Line &) const;

            bool operator <(const Point &) const;
            bool operator <=(const Point &) const;
            bool operator >(const Point &) const;
            bool operator >=(const Point &) const;

            friend bool operator <(const Point &, const Line &);
            friend bool operator <=(const Point &, const Line &);
            friend bool operator >(const Point &, const Line &);
            friend bool operator >=(const Point &, const Line &);

            // METHODS.

            double angular() const;
            bool contains(const Point &) const;
            bool is_parallel(const Line &) const;

            // OUTPUT.

            friend std::ostream &operator <<(std::ostream &, const Line &);
    };

    // Segment, implemented under src/Segment.cpp

    /**
     * @brief 2D Segment class.
     * 
     */
    class Segment {
        protected:

            // Extremes.
            const Point a;
            const Point b;

        public:

            // CONSTRUCTORS.

            Segment(const Point &, const Point &);
            Segment(const std::array<Point, 2> &);
            Segment(const Segment &);

            // READ.

            Point operator [](const std::size_t &) const;

            // COMPARISONS

            bool operator ==(const Segment &) const;

            // METHODS.

            Line line() const;

            bool contains(const Point &) const;
            bool contains(const Segment &) const;

            // OUTPUT.

            friend std::ostream &operator <<(std::ostream &, const Segment &);
    };

    // Polygon, implemented under src/Polygon.cpp

    /**
     * @brief 2D Polygon struct.
     * 
     */
    struct Polygon {

        // Points (Counterwise ordered).
        std::vector<Point> points;

        // CONSTRUCTORS.

        Polygon(const std::vector<Point> &);
        Polygon(const Polygon &);
        Polygon &operator =(const Polygon &);

        // METHODS.

        std::vector<Segment> edges() const;

        bool contains(const Point &) const;
        bool contains(const Segment &) const;

        double area() const;
        Point centroid() const;

        Point random() const;

        std::array<Point, 2> box() const;

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Polygon &);
    };

    // METHODS.
    // Implemented under src/Geometry.cpp

    Line bisector(const Point &, const Point &);
    
    std::vector<Point> intersections(const Line &, const Line &);
    std::vector<Point> intersections(const Line &, const Segment &);
    std::vector<Point> intersections(const Line &, const Polygon &);

    Polygon collapse(const Polygon &, const Point &);
    Polygon collapse(const Polygon &, const Segment &);

    Polygon reduce(const Polygon &, const Line &, const Point &);

}

namespace std {

    double abs(const pacs::Point &);
    double abs(const pacs::Segment &);

}

#endif