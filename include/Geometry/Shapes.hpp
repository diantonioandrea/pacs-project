/**
 * @file Shapes.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SHAPES_PACS
#define SHAPES_PACS

#include <Type.hpp>

#include <iostream>
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
            Real x;
            Real y;

        public:

            // CONSTRUCTORS.

            Point(const Real &, const Real &);
            Point(const std::array<Real, 2> &);
            Point(const Point &);
            Point &operator =(const Point &);

            // READ AND WRITE.

            Real operator [](const std::size_t &) const;
            Real &operator [](const std::size_t &);

            // COMPARISONS.

            bool operator ==(const Point &) const;
            bool operator !=(const Point &) const;

            // METHODS.

            bool is_zero() const;

            // OPERATORS.

            Point operator *(const Real &) const;
            friend Point operator *(const Real &, const Point &);
            Point &operator *=(const Real &);
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
            const Real a;
            const Real b;
            const Real c;

        public:

            // CONSTRUCTORS.

            Line(const Real &, const Real &, const Real &);
            Line(const std::array<Real, 3> &);
            Line(const Line &);

            // READ.

            Real operator [](const std::size_t &) const;
            Real x(const Real &) const;
            Real y(const Real &) const;

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

            bool contains(const Point &) const;

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

        Real area() const;
        Point centroid() const;

        Point random() const;

        std::array<Point, 2> box() const;

        // OUTPUT.

        friend std::ostream &operator <<(std::ostream &, const Polygon &);
    };

    // METHODS.
    // Implemented under src/Geometry.cpp

    // Distance between objects.
    Real distance(const Point &, const Point &);
    Real distance(const Point &, const Line &);
    Real distance(const Point &, const Segment &);

    Line bisector(const Point &, const Point &);
    
    // Intersections between objects.
    std::vector<Point> intersections(const Line &, const Line &);
    std::vector<Point> intersections(const Line &, const Segment &);
    std::vector<Point> intersections(const Line &, const Polygon &);

    // Normal line.
    Line normal(const Line &, const Point &);
    Line normal(const Segment &, const Point &);

    // Polygon collapse.
    Polygon collapse(const Polygon &, const Point &);
    Polygon collapse(const Polygon &, const Segment &);

    // Polygon Point reflection.
    std::vector<Point> reflections(const Polygon &, const Point &);

    // Polygon reduction.
    Polygon reduce(const Polygon &, const Line &, const Point &);

}

namespace std {

    pacs::Real abs(const pacs::Point &);
    pacs::Real abs(const pacs::Segment &);

}

#endif