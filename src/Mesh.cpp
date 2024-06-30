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
     * @brief Constructs a new Mesh from a given domain, diagram and degrees vector.
     * 
     * @param domain 
     * @param diagram 
     */
    Mesh::Mesh(const Polygon &domain, const std::vector<Polygon> &diagram, const std::vector<std::size_t> &degrees): domain{domain} {

        // Elements.
        this->elements = mesh_elements(diagram, degrees);

        // Neighbours.
        this->neighbours = mesh_neighbours(domain, this->elements);

        // Areas and biggest simplices.
        this->areas = mesh_areas(diagram);
        this->max_simplices = mesh_max_simplices(diagram);

        // Number of quadrature nodes and solution evaluation entries.
        std::size_t entries = 0;

        for(const auto &element: this->elements)
            entries += element.element.points.size();

        this->quadrature = 15; // Arbitrary.
        this->entries = entries * this->quadrature * this->quadrature;
    }

    /**
     * @brief Constructs a new Mesh from a given domain, diagram and uniform degree.
     * 
     * @param domain 
     * @param diagram 
     * @param degree 
     */
    Mesh::Mesh(const Polygon &domain, const std::vector<Polygon> &diagram, const std::size_t &degree): 
    Mesh(domain, diagram, std::vector<std::size_t>(diagram.size(), degree)) {}

    /**
     * @brief Copy constructor.
     * 
     * @param mesh 
     */
    Mesh::Mesh(const Mesh &mesh):
    domain{mesh.domain}, elements{mesh.elements}, neighbours{mesh.neighbours}, areas{mesh.areas}, max_simplices{std::vector<Vector<Real>>(mesh.max_simplices)}, quadrature{mesh.quadrature}, entries{mesh.entries} {}

    /**
     * @brief Copy operator.
     * 
     * @param mesh 
     * @return Mesh& 
     */
    Mesh &Mesh::operator =(const Mesh &mesh) {
        this->domain = mesh.domain;
        this->elements = mesh.elements;
        this->neighbours = mesh.neighbours;
        this->areas = mesh.areas;
        this->max_simplices = std::vector<Vector<Real>>(mesh.max_simplices);
        this->quadrature = mesh.quadrature;
        this->entries = mesh.entries;

        return *this;
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
    void Mesh::write(const std::string &filename, const bool &degrees) {
        // File loading.
        std::ofstream file{filename};

        file << "@ " << filename << "\n";
        file << "@ polyplot.py readable mesh\n";

        // Stats.
        file << "@ Domain: " << this->domain << "\n";
        file << "@ Elements number: " << this->elements.size() << "\n";

        // Polygons.
        file << "@ Elements: \n";

        for(const auto &element: this->elements) {
            Polygon polygon = element.element;

            for(const auto &vertex: polygon.points)
                file << std::setprecision(12) << vertex[0] << " " << std::setprecision(12) << vertex[1] << " ";
                        
            if(degrees) 
                file << element.degree;

            file << "\n";
        }
    }

}