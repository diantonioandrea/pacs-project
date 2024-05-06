/**
 * @file Mesh.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MESH_PACS
#define MESH_PACS

// Geometry.
#include <Geometry.hpp>

// Voronoi.
#include <Voronoi.hpp>

// Containers.
#include <vector>

namespace pacs {

    /**
     * @brief Element struct.
     * 
     */
    struct Element {
            
        // Indices.
        std::vector<std::size_t> nodes;
        std::vector<std::size_t> edges;

        // Polynomial degree.
        std::size_t degree;

        // CONSTRUCTORS.

        Element(const std::vector<std::size_t> &, const std::vector<std::size_t> &);

        Element(const std::vector<std::size_t> &, const std::vector<std::size_t> &, const std::size_t &);

    };

    /**
     * @brief Mesh class.
     * 
     */
    class Mesh {
        protected:

            // Geometry.
            std::vector<Point> nodes;
            std::vector<Segment> edges;

        public:

            // Elements.
            std::vector<Element> elements;

            // Indices.
            std::vector<std::size_t> boundary_nodes;
            std::vector<std::size_t> boundary_edges;

            // CONSTRUCTORS.

            Mesh(const Polygon &, const std::size_t &);
            Mesh(const Mesh &);

            // READ.
        
            Point node(const std::size_t &) const;
            Segment edge(const std::size_t &) const;
            Polygon element(const std::size_t &) const;
            Polygon element(const Element &) const;

            // STATS.

            std::size_t nodes_number() const;
            std::size_t edges_number() const;
            std::size_t elements_number() const;
    };

}

#endif