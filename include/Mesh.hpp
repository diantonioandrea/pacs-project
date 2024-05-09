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

// Type.
#include <Type.hpp>

// Geometry.
#include <Geometry.hpp>

// Voronoi.
#include <Voronoi.hpp>

// Vectors.
#include <Vector.hpp>

// Containers.
#include <vector>
#include <array>

// Output.
#include <string>

namespace pacs {

    // Element, implemented under Element.cpp

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

        // METHODS.

        inline std::size_t dofs() const { return (this->degree + 1) * (this->degree + 2) / 2; }

    };

    // Mesh, implemented under Mesh.cpp

    /**
     * @brief Mesh class.
     * 
     */
    class Mesh {
        protected:

            // Geometry.
            Polygon domain;
            std::vector<Point> nodes;
            std::vector<Segment> edges;

        public:

            // Elements.
            std::vector<Element> elements;

            // Indices.
            std::vector<std::size_t> boundary_nodes;
            std::vector<std::size_t> boundary_edges;

            // Neighbours.
            std::vector<std::vector<std::pair<std::size_t, int>>> neighbours;

            // Areas.
            std::vector<Real> areas;
            std::vector<Vector<Real>> max_simplices;

            // Penalty coefficient.
            Real penalty = 10.0;

            // CONSTRUCTORS.

            Mesh(const Polygon &, const std::size_t &);
            Mesh(const Polygon &, const std::vector<Polygon> &);

            // READ.
        
            Point node(const std::size_t &) const;
            Segment edge(const std::size_t &) const;
            Polygon element(const std::size_t &) const;
            Polygon element(const Element &) const;

            // STATS.

            inline std::size_t nodes_number() const { return this->nodes.size(); }
            inline std::size_t edges_number() const { return this->edges.size(); }
            inline std::size_t elements_number() const { return this->elements.size(); }
            std::size_t dofs() const;

            // OUTPUT.

            void write(const std::string &);
    };

    // METHODS.
    // Implemented under src/Builder.cpp

    std::vector<Polygon> mesh_diagram(const Polygon &, const std::size_t &);

    std::vector<Point> mesh_nodes(const std::vector<Polygon> &);
    std::vector<Segment> mesh_edges(const std::vector<Polygon> &);
    std::vector<Element> mesh_elements(const std::vector<Polygon> &, const std::vector<Point> &, const std::vector<Segment> &);

    std::vector<std::size_t> mesh_boundary_nodes(const Polygon &, const std::vector<Point> &);
    std::vector<std::size_t> mesh_boundary_edges(const Polygon &, const std::vector<Segment> &);

    std::vector<std::vector<std::pair<std::size_t, int>>> mesh_neighbours(const std::vector<Element> &, const std::vector<std::size_t> &);

    std::vector<Real> mesh_areas(const std::vector<Polygon> &);
    std::vector<Vector<Real>> mesh_max_simplices(const std::vector<Polygon> &);

}

#endif