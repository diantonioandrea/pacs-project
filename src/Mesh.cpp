/**
 * @file Mesh.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace pacs {

    // CONSTRUCTORS.

    /**
     * @brief Constructs a new Mesh from a given domain and diagram.
     * 
     * @param domain 
     * @param diagram 
     */
    Mesh::Mesh(const Polygon &domain, const std::vector<Polygon> &diagram, const std::size_t &degree): domain{domain} {

        // Elements.
        this->elements = mesh_elements(diagram, degree);

        // Neighbours.
        this->neighbours = mesh_neighbours(domain, this->elements);

        // Areas and biggest simplices.
        this->areas = mesh_areas(diagram);
        this->max_simplices = mesh_max_simplices(diagram);

        // Number of quadrature nodes and solution evaluation entries.
        std::size_t entries = 0;

        for(const auto &element: this->elements)
            entries += element.element.points.size();

        this->quadrature = 7; // Arbitrary.
        this->entries = entries * this->quadrature * this->quadrature;
    }

    // READ.

    /**
     * @brief Returns the j-th element's polygon.
     * 
     * @param j 
     * @return Polygon 
     */
    Polygon Mesh::element(const std::size_t &j) const {
        #ifndef NDEBUG // Integrity check.
        assert(j < this->elements.size());
        #endif

        return this->elements[j].element;
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
        file << "@ Elements: " << this->elements.size() << "\n";

        // Polygons.
        file << "@ Elements\' coordinates: \n";

        for(const auto &element: this->elements) {
            Polygon polygon = element.element;

            for(const auto &vertex: polygon.points)
                file << vertex[0] << " " << vertex[1] << " ";
            
            file << "\n";
        }
    }

}