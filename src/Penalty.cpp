/**
 * @file Penalty.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Fem.hpp>
#include <Geometry.hpp>

namespace pacs {

    /**
     * @brief Returns the penalty coefficients for a given element. "Max" policy.
     * 
     * @param mesh 
     * @param index 
     * @return Vector<Real> 
     */
    Vector<Real> penalty(const Mesh &mesh, const std::size_t &index) {
        
        // Element.
        Element element = mesh.elements[index];
        Polygon polygon = mesh.element(index);

        std::vector<std::array<int, 3>> neighbours = mesh.neighbours[index];

        // Sizes.
        std::vector<Segment> edges = polygon.edges();
        Vector<Real> sizes{element.edges.size()};

        for(std::size_t j = 0; j < sizes.length; ++j)
            sizes[j] = std::abs(edges[j]);

        // Element's area.
        Real area = mesh.areas[index];

        // Biggest simplices areas.
        Vector<Real> areas = mesh.max_simplices[index];

        // Inverse constant.
        Vector<Real> inverse = area / areas;

        // Coefficients.
        Real penalty_coefficient = mesh.penalty * (element.degree * element.degree);
        Vector<Real> penalty_dirichlet = penalty_coefficient * inverse * sizes / area;

        // Penalty evaluation.
        Vector<Real> penalties{neighbours.size()};
        Vector<Real> internal{neighbours.size()}; // Element.
        Vector<Real> external{neighbours.size()}; // Neighbour.
        Vector<Real> inverse_external{neighbours.size()};

        for(std::size_t j = 0; j < neighbours.size(); ++j) {
            if(neighbours[j][1] == -1) {
                penalties[j] = penalty_dirichlet[j];
                continue;
            }

            inverse_external[j] = mesh.areas[neighbours[j][1]] / mesh.max_simplices[neighbours[j][1]][neighbours[j][2]];
            internal[j] = penalty_coefficient * inverse[j] * sizes[j] / area;
            external[j] = penalty_coefficient * inverse_external[j] * sizes[j] / mesh.areas[neighbours[j][1]];

            penalties[j] = (internal[j] > external[j]) ? internal[j] : external[j];
        }

        return penalties;
    }

}