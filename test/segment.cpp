/**
 * @file segment.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Testing Segment.
#include <Geometry.hpp>

int main() {

    // Constructing some points.
    pacs::Point o{0.0, 0.0};
    pacs::Point a{1.0, 1.0};
    pacs::Point b{-1.0, 1.0};
    pacs::Point c{-1.0, -1.0};
    pacs::Point d{1.0, -1.0};

    // Constructing some Segments.
    pacs::Segment oa{o, a};
    pacs::Segment ob{o, b};
    pacs::Segment oc{o, c};
    pacs::Segment od{o, d};

    // Orientations output.
    std::cout << oa.orientation() << std::endl;
    std::cout << ob.orientation() << std::endl;
    std::cout << oc.orientation() << std::endl;
    std::cout << od.orientation() << std::endl;
    
}