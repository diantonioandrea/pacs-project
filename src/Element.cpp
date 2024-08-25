/**
 * @file Element.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <PacsHPDG.hpp>

#include <cassert>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Element from given nodes and edges.
     * 
     * @param nodes Element's vertices.
     * @param edges Element's edges.
     */
    Element::Element(const Polygon &element): element{element}, nodes{element.points}, edges{element.edges()}, degree{1} {}

    /**
     * @brief Construct a new Element:: Element object
     * 
     * @param nodes Element's vertices.
     * @param edges Element's edges
     * @param degree Element's degree.
     */
    Element::Element(const Polygon &element, const std::size_t &degree): element{element}, nodes{element.points}, edges{element.edges()}, degree{degree} {}

}