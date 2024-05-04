/**
 * @file line.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// Output.
#include<iostream>

// Testing Line.
#include <Geometry.hpp>

int main() {

    // Constructing two points.
    pacs::Point p{0.0, 1.0};
    pacs::Point q{1.0, 0.0};

    // Constructing two lines.
    pacs::Line bisector = pacs::bisector(p, q);
    pacs::Line line{1.0, 1.0, 1.0};

    // Intersection output.
    std::cout << pacs::intersections(line, bisector)[0] << std::endl;
    
}