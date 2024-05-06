/**
 * @file Mesh_Element.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Mesh.hpp>

// Assertions.
#include <cassert>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Element from given nodes and edges.
     * 
     * @param nodes 
     * @param edges 
     */
    Element::Element(const std::vector<std::size_t> &nodes, const std::vector<std::size_t> &edges): nodes{nodes}, edges{edges}, degree{1} {}

    /**
     * @brief Construct a new Element:: Element object
     * 
     * @param nodes 
     * @param edges 
     * @param degree 
     */
    Element::Element(const std::vector<std::size_t> &nodes, const std::vector<std::size_t> &edges, const std::size_t &degree): nodes{nodes}, edges{edges}, degree{degree} {}

    // METHODS.

    /**
     * @brief Returns the number of local degrees of freedom.
     * 
     * @return std::size_t 
     */
    inline std::size_t Element::dofs() const {
        return (this->degree + 1) * (this->degree + 2) / 2;
    }

}