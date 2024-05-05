/**
 * @file Mesh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Mesh.hpp>

// Assertions.
#include <cassert>

// Voronoi.
#include <Voronoi.hpp>

namespace pacs {

    // ELEMENT.

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Element from given nodes and edges.
     * 
     * @param nodes 
     * @param edges 
     */
    Element::Element(const std::vector<std::size_t> &nodes, const std::vector<std::size_t> &edges): nodes{nodes}, edges{edges} {}

    // MESH.

    // CONSTRUCTORS.

    /**
     * @brief Construct a new Mesh from a Voronoi diagram.
     * 
     * @param mesh 
     */
    Mesh::Mesh(const Polygon &domain, const std::size_t &cells) {
        
        // Diagram.
        std::vector<Polygon> mesh = voronoi(domain, cells);

        // Relaxation.
        for(std::size_t j = 0; j < 25; ++j)
            mesh = lloyd(domain, mesh);

        // Building nodes and edges.
        for(const auto &cell: mesh) {
            for(const auto &node: cell.vertices()) {
                bool flag = true;

                for(const auto &point: this->nodes) {
                    if(point == node) {
                        flag = false;
                        break;
                    }
                }

                if(flag)
                    this->nodes.emplace_back(node);
            }

            for(const auto &edge: cell.edges()) {
                bool flag = true;

                for(const auto &segment: this->edges) {
                    if(segment == edge) {
                        flag = false;
                        break;
                    }
                }

                if(flag)
                    this->edges.emplace_back(edge);
            }
        }

        // Elements.
        for(const auto &cell: mesh) {
            std::vector<std::size_t> nodes;
            std::vector<std::size_t> edges;

            for(const auto &node: cell.vertices()) {
                for(std::size_t j = 0; j < this->nodes.size(); ++j) {
                    if(node == this->nodes[j]) {
                        nodes.emplace_back(j);
                        break;
                    }
                }
            }

            for(const auto &edge: cell.edges()) {
                for(std::size_t j = 0; j < this->edges.size(); ++j) {
                    if(edge == this->edges[j]) {
                        edges.emplace_back(j);
                        break;
                    }
                }
            }

            this->elements.emplace_back(nodes, edges);
        }
        
        // Boundary nodes and edges.
        for(const auto &bound: domain.edges()) {
            for(std::size_t j = 0; j < this->nodes.size(); ++j) {
                if(bound.contains(this->nodes[j]))
                    this->boundary_nodes.emplace_back(j);
            }

            for(std::size_t j = 0; j < this->edges.size(); ++j) {
                if(bound.contains(this->edges[j]))
                    this->boundary_edges.emplace_back(j);
            }
        }

    }

    // READ.
        
    Point Mesh::node(const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->nodes.size());
        #endif

        return this->nodes[j];
    }

    Segment Mesh::edge(const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->edges.size());
        #endif

        return this->edges[j];
    }

    /**
     * @brief Copy constructor.
     * 
     * @param mesh 
     */
    Mesh::Mesh(const Mesh &mesh):
    nodes{mesh.nodes}, edges{mesh.edges}, elements{mesh.elements}, boundary_nodes{mesh.boundary_nodes}, boundary_edges{mesh.boundary_edges} {}

    // STATS.

    /**
     * @brief Returns the number of nodes.
     * 
     * @return std::size_t 
     */
    std::size_t Mesh::nodes_number() const { return this->nodes.size(); }

    /**
     * @brief Returns the number of edges.
     * 
     * @return std::size_t 
     */
    std::size_t Mesh::edges_number() const { return this->edges.size(); }

    /**
     * @brief Returns the number of elements.
     * 
     * @return std::size_t 
     */
    std::size_t Mesh::elements_number() const { return this->elements.size(); }

}