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
#include <array>

// Geomtry tolerance.
#ifndef GEOMETRY_TOLERANCE
#define GEOMETRY_TOLERANCE 1E-14
#endif

namespace pacs {
    
    /**
     * @brief 2D Point class.
     * 
     */
    class Point {
        private:

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

            // METHODS.

            double norm() const;
            bool is_origin() const;

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

    /**
     * @brief 2D Line (ax + by = c) class.
     * 
     */
    class Line {
        private:

            // Parameters.
            const double a;
            const double b;
            const double c;

        public:

            // CONSTRUCTORS.

            Line(const double &, const double &, const double &);
            Line(const std::array<double, 3> &);
            Line(const Line &);

            // OUTPUT.

            friend std::ostream &operator <<(std::ostream &, const Line &);
    };

    /**
     * @brief 2D Segment class.
     * 
     */
    class Segment {
        private:

            // Extremes.
            const Point a;
            const Point b;

        public:

            // CONSTRUCTORS.

            Segment(const Point &, const Point &);
            Segment(const std::array<Point, 2> &);
            Segment(const Segment &);

            // METHODS.

            Line line() const;

            // OUTPUT.

            friend std::ostream &operator <<(std::ostream &, const Segment &);
    };

    // METHODS.

    Line bisector(const Point &, const Point &);

}

#endif