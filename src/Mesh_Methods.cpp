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

// Small edges tolerance.
#ifndef COLLAPSE_TOLERANCE
#define COLLAPSE_TOLERANCE 1E-1
#endif

// Maximum number of iterations for the Lloyd's algorithm.
#ifndef LLOYD_MAX_ITER
#define LLOYD_MAX_ITER 256
#endif

namespace pacs {

    // METHODS.

    std::vector<Polygon> mesh_diagram(const Polygon &domain, const std::size_t &cells) {
        
        // Diagram.
        std::vector<Polygon> mesh = voronoi(domain, cells);

        // Relaxation.
        for(std::size_t j = 0; j < LLOYD_MAX_ITER; ++j)
            mesh = lloyd(domain, mesh);

        // Small edges collapse.
        double size = 0.0;

        for(std::size_t j = 0; j < mesh.size(); ++j) {
            for(const auto &edge: mesh[j].edges())
                size = (std::abs(edge[1] - edge[0]) > size) ? std::abs(edge[1] - edge[0]) : size;
        }

        std::size_t index = 0;
        while(index < mesh.size()) {
            bool flag = true;

            // Cannot collapse triangles.
            if(mesh[index].edges().size() == 3)
                continue;

            // Looks for small edges.
            for(auto &edge: mesh[index].edges()) {
                if(std::abs(edge[0] - edge[1]) > COLLAPSE_TOLERANCE * size)
                    continue;

                for(std::size_t k = 0; k < mesh.size(); ++k) {
                    if(index == k)
                        continue;

                    if(domain.contains(edge))
                        break;

                    // Collapses small edges.
                    if(mesh[k].contains(edge)) {
                        mesh[index] = collapse(mesh[index], edge);
                        mesh[k] = collapse(mesh[k], edge);
                        flag = false;

                        Point mid = (edge[0] + edge[1]) * 0.5;

                        // Moves vertices.
                        for(std::size_t h = 0; h < mesh.size(); ++h) {
                            if((index  == h) || (k == h))
                                continue;

                            for(std::size_t l = 0; l < mesh[h].points.size(); ++l) {
                                if((mesh[h].points[l] == edge[0]) || (mesh[h].points[l] == edge[1]))
                                    mesh[h].points[l] = mid;
                            }
                        }
                    }
                }

                if(!flag)
                    break;
            }

            if(flag)
                ++index;
        }

        return mesh;
    }

    /**
     * @brief Returns the vector of unique nodes inside a mesh.
     * 
     * @param mesh 
     * @return std::vector<Point> 
     */
    std::vector<Point> mesh_nodes(const std::vector<Polygon> &mesh) {
        std::vector<Point> nodes;

        for(const auto &cell: mesh) {
            for(const auto &candidate: cell.points) {
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

            for(const auto &node: cell.points) {
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

    std::vector<std::vector<std::pair<std::size_t, int>>> mesh_neighbours(const std::vector<Element> &elements, const std::vector<std::size_t> &boundary_edges) {
        std::vector<std::vector<std::pair<std::size_t, int>>> neighbours;

        for(std::size_t j = 0; j < elements.size(); ++j) {
            std::vector<std::pair<std::size_t, int>> element_neighbours;

            for(std::size_t k = 0; k < elements[j].edges.size(); ++k) {
                bool boundary = false;

                for(const auto &edge: boundary_edges) {
                    if(edge == elements[j].edges[k]) {
                        std::pair<std::size_t, int> pair{k, -1};
                        element_neighbours.emplace_back(pair);
                        boundary = true;
                        break;
                    }
                }

                if(boundary)
                    continue;

                for(std::size_t h = 0; h < elements.size(); ++h) {
                    if(j == h)
                        continue;

                    bool connection = false;

                    for(const auto &edge: elements[h].edges) {
                        if(edge == elements[j].edges[k]) {
                            std::pair<std::size_t, int> pair{k, h};
                            element_neighbours.emplace_back(pair);
                            connection = true;
                            break;
                        }
                    }

                    if(connection)
                        break;
                }
            }

            neighbours.emplace_back(element_neighbours);
        }

        return neighbours;
    }

}