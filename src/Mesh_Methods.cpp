/**
 * @file Mesh_Methods.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Mesh.hpp>

namespace pacs {

    // METHODS.

    /**
     * @brief Returns the vector of unique nodes inside a mesh.
     * 
     * @param mesh 
     * @return std::vector<Point> 
     */
    std::vector<Point> mesh_nodes(const std::vector<Polygon> &mesh) {
        std::vector<Point> nodes;

        for(const auto &cell: mesh) {
            for(const auto &candidate: cell.vertices()) {
                bool flag = true;

                for(const auto &node: nodes) {
                    if(candidate == node) {
                        flag = false;
                        break;
                    }
                }

                if(flag)
                    nodes.emplace_back(candidate);
            }
        }

        return nodes;
    }

    /**
     * @brief Returns the vector of unique edges inside a mesh.
     * 
     * @param mesh 
     * @return std::vector<Segment> 
     */
    std::vector<Segment> mesh_edges(const std::vector<Polygon> &mesh) {
        std::vector<Segment> edges;

        for(const auto &cell: mesh) {
            for(const auto &candidate: cell.edges()) {
                bool flag = true;

                for(const auto &edge: edges) {
                    if(candidate == edge) {
                        flag = false;
                        break;
                    }
                }

                if(flag)
                    edges.emplace_back(candidate);
            }
        }

        return edges;
    }
    
    /**
     * @brief Returns the vector of Elements inside a mesh.
     * 
     * @param mesh 
     * @param nodes 
     * @param edges 
     * @return std::vector<Element> 
     */
    std::vector<Element> mesh_elements(const std::vector<Polygon> &mesh, const std::vector<Point> &nodes, const std::vector<Segment> &edges) {
        std::vector<Element> elements;

        for(const auto &cell: mesh) {
            std::vector<std::size_t> element_nodes;
            std::vector<std::size_t> element_edges;

            for(const auto &node: cell.vertices()) {
                for(std::size_t j = 0; j < nodes.size(); ++j) {
                    if(node == nodes[j]) {
                        element_nodes.emplace_back(j);
                        break;
                    }
                }
            }

            for(const auto &edge: cell.edges()) {
                for(std::size_t j = 0; j < edges.size(); ++j) {
                    if(edge == edges[j]) {
                        element_edges.emplace_back(j);
                        break;
                    }
                }
            }

            elements.emplace_back(element_nodes, element_edges);
        }

        return elements;
    }

    /**
     * @brief Returns the vector of boundary nodes in a mesh given the domain.
     * 
     * @param domain 
     * @param nodes 
     * @return std::vector<std::size_t> 
     */
    std::vector<std::size_t> mesh_boundary_nodes(const Polygon &domain, const std::vector<Point> &nodes) {
        std::vector<std::size_t> boundary_nodes;

        for(const auto &bound: domain.edges()) {
            for(std::size_t j = 0; j < nodes.size(); ++j) {
                if(bound.contains(nodes[j]))
                    boundary_nodes.emplace_back(j);
            }
        }

        return boundary_nodes;
    }

    /**
     * @brief Returns the vector of boundary edges in a mesh given the domain.
     * 
     * @param domain 
     * @param edges 
     * @return std::vector<std::size_t> 
     */
    std::vector<std::size_t> mesh_boundary_edges(const Polygon &domain, const std::vector<Segment> &edges) {
        std::vector<std::size_t> boundary_edges;

        for(const auto &bound: domain.edges()) {
            for(std::size_t j = 0; j < edges.size(); ++j) {
                if(bound.contains(edges[j]))
                    boundary_edges.emplace_back(j);
            }
        }

        return boundary_edges;
    }

}