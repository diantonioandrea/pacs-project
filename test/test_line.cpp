/**
 * @file line.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Pacs.hpp>

#include <iostream>

int main() {

    // Constructing two points.
    pacs::Point p{0.0, 1.0};
    pacs::Point q{1.0, 0.0};

    // Constructing some Lines.
    pacs::Line bisector = pacs::bisector(p, q);
    pacs::Line bisector_symmetric = pacs::bisector(q, p);
    pacs::Line line{1.0, 1.0, 1.0};

    // Bisector output test.
    std::cout << bisector << std::endl;
    std::cout << bisector_symmetric << std::endl;

    // Intersection output.
    std::cout << pacs::intersections(line, bisector)[0] << std::endl;
    
}