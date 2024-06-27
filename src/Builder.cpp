/**
 * @file Builder.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-05-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Geometry.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace pacs {

    // METHODS.

    /**
     * @brief Returns a mesh diagram of a given domain.
     * 
     * @param domain 
     * @param cells 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> mesh_diagram(const Polygon &domain, const std::size_t &cells, const bool &uniform, const bool &reflect) {
        
        #ifndef NVERBOSE
        std::cout << "Generating a diagram for: " << domain << std::endl;
        #endif

        // Diagram.
        std::vector<Polygon> diagram;

        if(uniform)
            diagram = voronoi_uniform(domain, cells, reflect);
        else
            diagram = voronoi_random(domain, cells, reflect);

        #ifndef NVERBOSE
        std::cout << "\tGenerated the Voronoi diagram." << std::endl;
        #endif

        // Relaxation.
        diagram = mesh_relax(domain, diagram, reflect);

        // Reflection fixes.
        if(reflect) {
            for(const auto &vertex: domain.points) {
                for(std::size_t j = 0; j < diagram.size(); ++j) {
                    std::vector<Segment> edges_j = diagram[j].edges();

                    for(std::size_t k = 0; k < diagram[j].points.size(); ++k) {
                        if((edges_j[k].contains(vertex)) && (edges_j[k][0] != vertex) && (edges_j[k][1] != vertex)) {

                            Segment previous = edges_j[(k > 0) ? k - 1 : diagram[j].points.size() - 1];
                            Segment next = edges_j[(k < diagram[j].points.size() - 1) ? k + 1 : 0];

                            for(std::size_t h = 0; h < diagram.size(); ++h) {
                                if(h == j)
                                    continue;

                                std::vector<Segment> edges_h = diagram[h].edges();

                                for(std::size_t l = 0; l < edges_h.size(); ++l) {
                                    if(edges_h[l].contains(previous)) {
                                        if((edges_h[l][0] != previous[0]) && (edges_h[l][0] != previous[1]))
                                            diagram[h].points[l] = vertex;

                                        if((edges_h[l][1] != previous[0]) && (edges_h[l][1] != previous[1]))
                                            diagram[h].points[(l < edges_h.size() - 1) ? l + 1 : 0] = vertex;

                                        break;
                                    } else if(edges_h[l].contains(next)) {
                                        if((edges_h[l][0] != next[0]) && (edges_h[l][0] != next[1]))
                                            diagram[h].points[l] = vertex;

                                        if((edges_h[l][1] != next[0]) && (edges_h[l][1] != next[1]))
                                            diagram[h].points[(l < edges_h.size() - 1) ? l + 1 : 0] = vertex;

                                        break;
                                    }
                                }
                            }

                            diagram[j].points.erase(diagram[j].points.begin() + k);
                            diagram[j].points[k] = vertex;
                            break;

                        }
                    }
                }
            }
        }

        // Small edges collapse.
        std::vector<Real> sizes;
        sizes.resize(diagram.size(), 0.0);

        for(std::size_t j = 0; j < diagram.size(); ++j)
            for(const auto &edge: diagram[j].edges())
                sizes[j] = (std::abs(edge[1] - edge[0]) > sizes[j]) ? std::abs(edge[1] - edge[0]) : sizes[j];

        #ifndef NVERBOSE
        std::cout << "Evaluated elements sizes." << std::endl;
        #endif

        std::size_t index = 0;
        while(index < diagram.size()) {
            bool flag = true;

            #ifndef NVERBOSE
            std::cout << "\tCollapsing element " << index << "." << std::endl;
            #endif

            // Cannot collapse triangles.
            if(diagram[index].edges().size() == 3)
                continue;

            // Looks for small edges.
            for(auto &edge: diagram[index].edges()) {
                if(std::abs(edge) > COLLAPSE_TOLERANCE * sizes[index])
                    continue;

                if(domain.contains(edge))
                    continue;

                for(std::size_t k = 0; k < diagram.size(); ++k) {
                    if(index == k)
                        continue;

                    // Collapses small edges.
                    if(diagram[k].contains(edge)) {
                        diagram[index] = collapse(diagram[index], edge);
                        diagram[k] = collapse(diagram[k], edge);
                        flag = false;

                        // Collapse target.
                        Point target = (edge[0] + edge[1]) * 0.5;

                        // Avoids domain distortion.
                        for(const auto &boundary_edge: domain.edges()) {
                            if(boundary_edge.contains(edge[0])) {
                                target = edge[0];
                                break;
                            }

                            if(boundary_edge.contains(edge[1])) {
                                target = edge[1];
                                break;
                            }
                        }

                        // Moves vertices.
                        for(std::size_t h = 0; h < diagram.size(); ++h) {
                            if((index  == h) || (k == h))
                                continue;

                            for(std::size_t l = 0; l < diagram[h].points.size(); ++l) {
                                if((diagram[h].points[l] == edge[0]) || (diagram[h].points[l] == edge[1]))
                                    diagram[h].points[l] = target;
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

        return diagram;
    }

    /**
     * @brief Reads a diagram from a file.
     * 
     * @param filename 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> mesh_diagram(const std::string &filename) {

        #ifndef NVERBOSE
        std::cout << "Loading a diagram from: " << filename << std::endl;
        #endif

        // File loading.
        std::ifstream file{filename};

        // Diagram.
        std::vector<Polygon> diagram;

        // Reading.
        std::string line;
        
        while (std::getline(file, line)) {

            // Skip lines starting with '@'
            if (!line.empty() && line[0] == '@') {
                continue;
            }

            // Reading points.
            std::istringstream lineStream{line};

            std::vector<Point> points;
            Real x, y;

            while (lineStream >> x >> y) {
                Point point{x, y};
                points.emplace_back(point);
            }

            diagram.emplace_back(points);
        }

        return diagram;
    }

    /**
     * @brief Relaxes a diagram through Lloyd's algorithm.
     * 
     * @param domain 
     * @param elements 
     * @return std::vector<Polygon> 
     */
    std::vector<Polygon> mesh_relax(const Polygon &domain, const std::vector<Polygon> &elements, const bool &reflect) {

        std::vector<Polygon> diagram{elements};

        // Relaxation, Lloyd's algorithm.
        std::vector<Point> centroids(diagram.size(), Point{0.0, 0.0});

        for(std::size_t j = 0; j < LLOYD_MAX_ITER; ++j) {

            // Update.
            Real residual = 0.0;

            for(std::size_t k = 0; k < diagram.size(); ++k) {
                
                // New point.
                Point centroid = diagram[k].centroid();
                Real shift = distance(centroid, centroids[k]);

                // Residual (biggest centroid shift).
                residual = (shift > residual) ? shift : residual;

                // Update.
                centroids[k] = centroid;

            }

            // Relaxation step.
            diagram = voronoi(domain, centroids, reflect);

            #ifndef NVERBOSE
            std::cout << "\tCompleted step " << j + 1 << " of the Lloyd's algorithm. Residual: " << residual << std::endl;
            #endif

            if(residual < LLOYD_TOLERANCE)
                break;
        }

        return diagram;
    }

    /**
     * @brief Refines specified elements' size for a mesh.
     * 
     * @param mesh 
     * @param mask 
     */
    void mesh_refine_size(Mesh &mesh, const Mask &mask) {
        #ifndef NDEBUG // Integrity check.
        assert(mask.size() == mesh.elements.size());
        #endif

        #ifndef NVERBOSE
        std::cout << "Refining mesh size." << std::endl;
        #endif

        // Degrees.
        std::vector<std::size_t> degrees;

        for(std::size_t j = 0; j < mask.size(); ++j)
            if(!(mask[j]))
                degrees.emplace_back(mesh.elements[j].degree);

        // Polygons.
        std::vector<Polygon> diagram;
        std::vector<Polygon> refine;
        std::vector<Polygon> refined;

        for(std::size_t j = 0; j < mask.size(); ++j)
            if(mask[j]) {
                refine.emplace_back(mesh.element(j));

                for(std::size_t k = 0; k < mesh.elements[j].edges.size() + static_cast<std::size_t>(mesh.elements[j].edges.size() > 4); ++k)
                    degrees.emplace_back(mesh.elements[j].degree);

            } else
                diagram.emplace_back(mesh.element(j));

        // Refine.
        for(const auto &polygon: refine) {
            Point centroid = polygon.centroid();
            std::vector<Point> points;

            for(const auto &edge: polygon.edges()) {

                // New point.
                Point point = (edge[0] + edge[1]) * 0.5;
                points.emplace_back(point);

                for(auto &element: diagram) {
                    std::vector<Segment> edges = element.edges();

                    for(std::size_t k = 0; k < edges.size(); ++k)
                        if(edges[k] == edge) {

                            // Diagram editing.
                            element.points.insert(element.points.begin() + k + 1, point);
                            break;
                        }
                }
            }

            // New polygons.
            std::vector<Segment> edges = polygon.edges();
            std::vector<Point> centrals;

            for(std::size_t j = 0; j < points.size(); ++j) {
                Point first = points[j];
                Point second = (j < points.size() - 1) ? points[j + 1] : points[0];

                std::vector<Point> vertices{first, edges[j][1], second};

                if(edges.size() <= 4)
                    vertices.emplace_back(centroid);
                else {
                    vertices.emplace_back((second + centroid) * 0.5);
                    vertices.emplace_back((first + centroid) * 0.5);
                    centrals.emplace_back((second + centroid) * 0.5);
                }

                refined.emplace_back(Polygon{vertices});
            }

            // Central Polygon.
            if(edges.size() > 4)
                refined.emplace_back(Polygon{centrals});
        }

        // Update.
        for(const auto &polygon: refined)
            diagram.emplace_back(polygon);

        // Refinement.
        mesh = Mesh{mesh.domain, diagram, degrees};
    }

    /**
     * @brief Refines specified elements' degree for a mesh.
     * 
     * @param mesh 
     * @param mask 
     */
    void mesh_refine_degree(Mesh &mesh, const Mask &mask) {
        #ifndef NDEBUG // Integrity check.
        assert(mask.size() == mesh.elements.size());
        #endif

        #ifndef NVERBOSE
        std::cout << "Refining mesh degree." << std::endl;
        #endif

        for(std::size_t j = 0; j < mask.size(); ++j)
            if(mask[j])
                if(mesh.elements[j].degree < 6) // Artificial limitation. [!]
                    ++mesh.elements[j].degree;
    }
    
    /**
     * @brief Returns the vector of Elements inside a mesh.
     * 
     * @param mesh 
     * @param nodes 
     * @param edges 
     * @return std::vector<Element> 
     */
    std::vector<Element> mesh_elements(const std::vector<Polygon> &diagram, const std::vector<std::size_t> &degrees) {
        #ifndef NDEBUG
        assert(degrees.size() == diagram.size());
        #endif

        std::vector<Element> elements;

        #ifndef NVERBOSE
        std::cout << "Evaluating mesh elements." << std::endl;
        #endif

        for(std::size_t j = 0; j < diagram.size(); ++j)
            elements.emplace_back(diagram[j], degrees[j]);

        return elements;
    }

    /**
     * @brief Returns the neighbouring structure.
     * 
     * @param domain 
     * @param elements 
     * @return std::vector<std::vector<std::array<int, 3>>> 
     */
    std::vector<std::vector<std::array<int, 3>>> mesh_neighbours(const Polygon &domain, const std::vector<Element> &elements) {
        std::vector<std::vector<std::array<int, 3>>> neighbours;
        neighbours.resize(elements.size());

        #ifndef NVERBOSE
        std::cout << "Evaluating mesh neighbours." << std::endl;
        #endif

        for(std::size_t j = 0; j < elements.size(); ++j) {
            Element element = elements[j];

            std::vector<std::array<int, 3>> element_neighbours;
            element_neighbours.resize(element.edges.size());

            for(std::size_t k = 0; k < element.edges.size(); ++k) {

                // Connection.
                bool connection = false;

                for(std::size_t h = 0; h < elements.size(); ++h) {
                    if(j == h)
                        continue;

                    Element other = elements[h];

                    for(std::size_t l = 0; l < other.edges.size(); ++l)
                        if(element.edges[k] == other.edges[l]) {
                            element_neighbours[k][0] = static_cast<int>(k);
                            element_neighbours[k][1] = static_cast<int>(h);
                            element_neighbours[k][2] = static_cast<int>(l);
                            connection = true;
                            break;
                        }

                    if(connection)
                        break;
                }

                // Boundary.
                if(!connection) {
                    element_neighbours[k][0] = static_cast<int>(k);
                    element_neighbours[k][1] = -1;
                    element_neighbours[k][2] = -1;
                }
            }

            neighbours[j] = element_neighbours;
        }

        return neighbours;
    }

    /**
     * @brief Returns a vector of the elements' areas.
     * 
     * @param polygons 
     * @return std::vector<Real> 
     */
    std::vector<Real> mesh_areas(const std::vector<Polygon> &polygons) {
        std::vector<Real> areas;

        #ifndef NVERBOSE
        std::cout << "Evaluating elements' areas." << std::endl;
        #endif

        for(std::size_t j = 0; j < polygons.size(); ++j)
            areas.emplace_back(polygons[j].area());

        return areas;
    }

    /**
     * @brief Returns a vector of the biggest simplices area for every element's edge.
     * 
     * @param polygons 
     * @return std::vector<Vector<Real>> 
     */
    std::vector<Vector<Real>> mesh_max_simplices(const std::vector<Polygon> &polygons) {
        std::vector<Vector<Real>> max_simplices;

        #ifndef NVERBOSE
        std::cout << "Evaluating elements' biggest simplices." << std::endl;
        #endif

        for(const auto &polygon: polygons) {
            Vector<Real> areas{polygon.edges().size()};

            for(std::size_t j = 0; j < areas.length; ++j) {
                Segment edge = polygon.edges()[j];

                for(std::size_t k = 0; k < polygon.points.size(); ++k) {
                    if((polygon.points[k] == edge[0]) || (polygon.points[k] == edge[1]))
                        continue;

                    Polygon triangle{{edge[0], edge[1], polygon.points[k]}};
                    Real area = triangle.area();

                    areas[j] = (area > areas[j]) ? area : areas[j];
                }
            }

            max_simplices.emplace_back(areas);
        }

        return max_simplices;
    }

}