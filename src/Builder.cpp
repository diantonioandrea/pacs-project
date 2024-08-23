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
    std::vector<Polygon> mesh_diagram(const Polygon &domain, const std::size_t &cells, const bool &reflect, const bool &uniform) {
        
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

                // On-edge vertex.
                for(std::size_t j = 0; j < diagram.size(); ++j) {
                    std::vector<Segment> edges_j = diagram[j].edges();

                    for(std::size_t ej = 0; ej < edges_j.size(); ++ej) {

                        // Edge.
                        Segment edge_j = edges_j[ej];

                        if(!(edge_j.contains(vertex)) || (edge_j[0] == vertex) || (edge_j[1] == vertex))
                            continue;

                        // Pathologies.
                        std::size_t counter = 0;
                        std::size_t index_j = ej;
                        std::vector<std::array<std::size_t, 2>> edges;

                        for(std::size_t k = 0; k < diagram.size(); ++k) {
                            if(j == k)
                                continue;

                            std::vector<Segment> edges_k = diagram[k].edges();

                            for(std::size_t ek = 0; ek < edges_k.size(); ++ek) {

                                // Edge.
                                Segment edge_k = edges_k[ek];

                                if((edge_k.contains(edge_j[0])) && (edge_k[0] != edge_j[0]) && (edge_k[1] != edge_j[0])) {
                                    edges.emplace_back(std::array<std::size_t, 2>{k, ek});
                                    ++counter;
                                    break;
                                }

                                if((edge_k.contains(edge_j[1])) && (edge_k[0] != edge_j[1]) && (edge_k[1] != edge_j[1])) {
                                    edges.emplace_back(std::array<std::size_t, 2>{k, ek});
                                    index_j = (ej == diagram[j].points.size() - 1) ? 0 : ej + 1;
                                    ++counter;
                                    break;
                                }
                            }
                        }

                        // Controlled pathology.
                        assert(counter <= 2);
                        
                        // Neighbours' fix.
                        for(auto &[k, ek]: edges) {
                            std::size_t index_k_0 = ek;
                            std::size_t index_k_1 = (ek == diagram[k].points.size() - 1) ? 0 : ek + 1;

                            // Move.
                            bool moved = false;

                            for(const auto &edge: edges_j)
                                if(edge.contains(diagram[k].points[index_k_0])) {
                                    diagram[k].points[index_k_1] = vertex;
                                    moved = true;
                                    break;
                                }

                            if(moved)
                                continue;

                            for(const auto &edge: edges_j)
                                if(edge.contains(diagram[k].points[index_k_1]))
                                    diagram[k].points[index_k_0] = vertex;
                        }

                        if(counter == 1) { // Single error.

                            diagram[j].points[index_j] = vertex;

                        } else if(counter == 2) { // Double error.

                            diagram[j].points.erase(diagram[j].points.begin() + ej);
                            diagram[j].points[ej] = vertex;
                        }
                    }
                }
            }
        }

        // Small edges collapse.
        std::vector<Real> sizes;
        sizes.resize(diagram.size(), 0.0);

        for(std::size_t j = 0; j < diagram.size(); ++j)
            for(const auto &p: diagram[j].points)
                for(const auto &q: diagram[j].points)
                    sizes[j] = (distance(p, q) > sizes[j]) ? distance(p, q) : sizes[j];

        #ifndef NVERBOSE
        std::cout << "Collapsing small edges." << std::endl;
        #endif

        for(std::size_t j = 0; j < diagram.size(); ++j) {
            std::vector<Segment> edges_j = diagram[j].edges();

            for(std::size_t ej = 0; ej < edges_j.size(); ++ej) {

                // Edge.
                Segment edge_j = edges_j[ej];

                // Triangles.
                if(edges_j.size() == 3)
                    continue;

                // Boundary edge.
                if(domain.contains(edge_j))
                    continue;

                // Long edge.
                if(std::abs(edge_j) > COLLAPSE_TOLERANCE * sizes[j])
                    continue;

                // "Boundary" edge.
                bool bound = false;

                for(const auto &edge_d: domain.edges())
                    if((edge_d.contains(edge_j[0])) || (edge_d.contains(edge_j[1])))
                        bound = true;

                if(bound)
                    continue;

                #ifndef NVERBOSE
                std::cout << "\tCollapsing element " << j << "." << std::endl;
                #endif
                
                // Midpoint.
                Point midpoint = (edge_j[0] + edge_j[1]) * 0.5;

                // Other elements fix.
                for(std::size_t k = 0; k < diagram.size(); ++k) {
                    if(j == k)
                        continue;

                    std::vector<Segment> edges_k = diagram[k].edges();

                    // Same edge.
                    bool edge = false;

                    for(std::size_t ek = 0; ek < edges_k.size(); ++ek) {

                        // Edge.
                        Segment edge_k = edges_k[ek];

                        if(edge_j == edge_k) {
                            edge = true;

                            // Edit.
                            if(ek == edges_k.size() - 1) {
                                diagram[k].points.erase(diagram[k].points.begin() + edges_k.size() - 1);
                                diagram[k].points[0] = midpoint;
                            } else {
                                diagram[k].points.erase(diagram[k].points.begin() + ek + 1);
                                diagram[k].points[ek] = midpoint;
                            }

                            break;
                        }
                    }

                    if(edge)
                        continue;

                    // Same point.
                    for(std::size_t pk = 0; pk < diagram[k].points.size(); ++pk)
                        if((diagram[k].points[pk] == edge_j[0]) || (diagram[k].points[pk] == edge_j[1])) {
                            diagram[k].points[pk] = midpoint;
                            break;
                        }
                }

                // Element fix.
                if(ej == edges_j.size() - 1) {
                    diagram[j].points.erase(diagram[j].points.begin() + edges_j.size() - 1);
                    diagram[j].points[0] = midpoint;
                } else {
                    diagram[j].points.erase(diagram[j].points.begin() + ej + 1);
                    diagram[j].points[ej] = midpoint;
                }

                break;
            }
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

        for(std::size_t j = 0; j < LLOYD_ITER_MAX; ++j) {

            // Update.
            Real residual = 0.0;

            for(std::size_t k = 0; k < diagram.size(); ++k) {
                
                // New point.
                Point centroid = diagram[k].centroid();

                // Residual (Centroid shift).
                residual += distance(centroid, centroids[k]);

                // Update.
                centroids[k] = centroid;

            }

            // Relaxation step.
            diagram = voronoi(domain, centroids, reflect);

            #ifndef NVERBOSE
            std::cout << "\tCompleted step " << j + 1 << " of the Lloyd's algorithm. Residual: " << residual << std::endl;
            #endif

            if(residual < LLOYD_TOLERANCE * diagram.size())
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
        std::cout << "Refining mesh size: " << std::flush;
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

        #ifndef NVERBOSE
        std::cout << refine.size() << " elements." << std::endl;
        #endif

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
        std::cout << "Refining mesh degree: " << std::flush;
        std::size_t counter = 0;
        #endif

        for(std::size_t j = 0; j < mask.size(); ++j)
            if(mask[j]) {
                ++mesh.elements[j].degree;
                
                #ifndef NVERBOSE
                ++counter;
                #endif
            }

        #ifndef NVERBOSE
        std::cout << counter << " elements." << std::endl;
        #endif
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