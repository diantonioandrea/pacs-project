/**
 * @file test_fit.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-07-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <iostream>

int main() {

    // Constructing some data.
    pacs::Vector<pacs::Real> x{4};
    pacs::Vector<pacs::Real> y{4};

    x[0] = 1.0L;
    x[1] = 2.0L;
    x[2] = 3.0L;
    x[3] = 4.0L;

    y[0] = 1.0L;
    y[1] = 4.0L;
    y[2] = 9.0L;
    y[3] = 16.0L;

    // Polynomial fit.
    std::cout << pacs::polyfit(x, y, 2) << std::endl;
}