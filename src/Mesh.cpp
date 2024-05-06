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

// OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// IO handling.
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

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
    
    /**
     * @brief Returns the j-th node.
     * 
     * @param j 
     * @return Point 
     */
    Point Mesh::node(const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->nodes.size());
        #endif

        return this->nodes[j];
    }

    /**
     * @brief Returns the j-th edge.
     * 
     * @param j 
     * @return Segment 
     */
    Segment Mesh::edge(const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->edges.size());
        #endif

        return this->edges[j];
    }

    /**
     * @brief Returns the j-th element.
     * 
     * @param j 
     * @return Polygon 
     */
    Polygon Mesh::element(const std::size_t &j) const {
        std::vector<Point> nodes;

        for(const auto &index: this->elements[j].nodes)
            nodes.emplace_back(this->nodes[index]);

        return Polygon{nodes};
    }

    /**
     * @brief Returns an element.
     * 
     * @param element 
     * @return Polygon 
     */
    Polygon Mesh::element(const Element &element) const {
        std::vector<Point> nodes;

        for(const auto &index: element.nodes)
            nodes.emplace_back(this->nodes[index]);

        return Polygon{nodes};
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
    inline std::size_t Mesh::nodes_number() const { return this->nodes.size(); }

    /**
     * @brief Returns the number of edges.
     * 
     * @return std::size_t 
     */
    inline std::size_t Mesh::edges_number() const { return this->edges.size(); }

    /**
     * @brief Returns the number of elements.
     * 
     * @return std::size_t 
     */
    inline std::size_t Mesh::elements_number() const { return this->elements.size(); }

    /**
     * @brief Returns the number of degrees of freedom.
     * 
     * @return std::size_t 
     */
    inline std::size_t Mesh::dofs() const {
        std::size_t dofs = 0;

        #pragma omp parallel for reduction(+: dofs)
        for(const auto &element: this->elements)
            dofs += element.dofs();

        return dofs;
    }

    // OUTPUT.

    /**
     * @brief Outputs the mesh to a polyplot.py readable file.
     * 
     * @param filename 
     */
    void Mesh::write(const std::string &filename) {
        // File loading.
        std::ofstream file{filename};

        file << "@ polyplot.py readable mesh\n";

        // Stats.
        file << "@ Elements: " << this->elements_number() << "\n";
        file << "@ Nodes: " << this->elements_number() << "\n";
        file << "@ Edges: " << this->elements_number() << "\n";

        // Polygons.
        file << "@ Elements\' coordinates: \n";

        for(const auto &element: this->elements)
            file << this->element(element) << "\n";
    }

}