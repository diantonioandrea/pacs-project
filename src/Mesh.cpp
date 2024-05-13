/**
 * @file Mesh.cpp
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

// OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// IO handling.
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

// Voronoi.
#include <Voronoi.hpp>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Mesh from a given domain and diagram.
     * 
     * @param domain 
     * @param mesh 
     */
    Mesh::Mesh(const Polygon &domain, const std::vector<Polygon> &mesh): domain{domain} {

        // Building nodes and edges.
        this->nodes = mesh_nodes(mesh);
        this->edges = mesh_edges(mesh);

        // Elements.
        this->elements = mesh_elements(mesh, this->nodes, this->edges);
        
        // Boundary nodes and edges.
        this->boundary_nodes = mesh_boundary_nodes(domain, this->nodes);
        this->boundary_edges = mesh_boundary_edges(domain, this->edges);

        // Neighbours.
        this->neighbours = mesh_neighbours(this->elements, this->boundary_edges);

        // Areas and biggest simplices.
        this->areas = mesh_areas(mesh);
        this->max_simplices = mesh_max_simplices(mesh);

        // Degree and entries.
        std::size_t degree = 0, entries = 0;

        for(const auto &element: this->elements) {
            degree = (degree < element.degree) ? element.degree : degree;
            entries += element.edges.size();
        }

        entries *= degree * degree;

        this->degree = degree;
        this->entries = entries;
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

    // STATS.

    /**
     * @brief Returns the number of degrees of freedom.
     * 
     * @return std::size_t 
     */
    std::size_t Mesh::dofs() const {
        std::size_t dofs = 0;

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

        file << "@ " << filename << "\n";
        file << "@ polyplot.py readable mesh\n";

        // Stats.
        file << "@ Domain: " << this->domain << "\n";
        file << "@ Elements: " << this->elements_number() << "\n";
        file << "@ Nodes: " << this->nodes_number() << "\n";
        file << "@ Edges: " << this->edges_number() << "\n";

        // Polygons.
        file << "@ Elements\' coordinates: \n";

        for(const auto &element: this->elements)
            file << this->element(element) << "\n";
    }

}