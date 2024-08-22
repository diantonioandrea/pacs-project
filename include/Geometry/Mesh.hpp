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

#include <Base.hpp>
#include <Algebra.hpp>

#include "Shapes.hpp"
#include "Shapes.hpp"

#include <vector>
#include <array>
#include <string>

#ifndef COLLAPSE_TOLERANCE
#define COLLAPSE_TOLERANCE 1E-1
#endif

#ifndef LLOYD_TOLERANCE
#define LLOYD_TOLERANCE 1E-4
#endif

#ifndef LLOYD_MAX_ITER
#define LLOYD_MAX_ITER 2E2
#endif

namespace pacs {

    // Element, implemented under Element.cpp

    /**
     * @brief Element struct. Polygon wrapper.
     * 
     */
    struct Element {
            
        // Polygon.
        Polygon element;
        
        std::vector<Point> nodes;
        std::vector<Segment> edges;

        // Polynomial degree.
        std::size_t degree;

        // CONSTRUCTORS.

        Element(const Polygon &);
        Element(const Polygon &, const std::size_t &);

        // METHODS.

        inline std::size_t dofs() const { return (this->degree + 1) * (this->degree + 2) / 2; }

    };

    // Mesh, implemented under Mesh.cpp

    /**
     * @brief Mesh struct.
     * 
     */
    struct Mesh {

        // Geometry.
        Polygon domain; // Polygonal domain.

        // Elements.
        std::vector<Element> elements;

        // Neighbours.
        std::vector<std::vector<std::array<int, 3>>> neighbours;

        // Areas.
        std::vector<Real> areas;
        std::vector<Vector<Real>> max_simplices;

        // Penalty coefficient.
        Real penalty = 10.0;

        // Quadrature nodes.
        std::size_t quadrature;

        // Entries for the solution.
        std::size_t entries;

        // CONSTRUCTORS.

        Mesh(const Polygon &, const std::vector<Polygon> &, const std::vector<std::size_t> &degrees);
        Mesh(const Polygon &, const std::vector<Polygon> &, const std::size_t &degree = 1);
        Mesh(const Mesh &);

        Mesh &operator =(const Mesh &);

        // READ, WRAPPERS.
    
        Polygon element(const std::size_t &) const;

        // STATS.

        std::size_t dofs() const;

        // OUTPUT.

        void write(const std::string &, const bool &degrees = false);
    };

    // METHODS.
    // Implemented under src/Builder.cpp

    // Diagrams.
    std::vector<Polygon> mesh_diagram(const Polygon &, const std::size_t &, const bool &uniform = false, const bool &reflect = false);
    std::vector<Polygon> mesh_diagram(const std::string &);
    std::vector<Polygon> mesh_relax(const Polygon &, const std::vector<Polygon> &, const bool &reflect = false);

    // Refinement.
    void mesh_refine_size(Mesh &, const Mask &);
    void mesh_refine_degree(Mesh &, const Mask &);

    // Data.
    std::vector<Element> mesh_elements(const std::vector<Polygon> &, const std::vector<std::size_t> &);
    std::vector<std::vector<std::array<int, 3>>> mesh_neighbours(const Polygon &, const std::vector<Element> &);

    std::vector<Real> mesh_areas(const std::vector<Polygon> &);
    std::vector<Vector<Real>> mesh_max_simplices(const std::vector<Polygon> &);

}

#endif