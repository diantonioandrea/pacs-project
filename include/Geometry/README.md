# `include/Geometry/`

## Introduction

This directory contains the implementation of geometry-related components for the project.

## Classes and structs

### [`include/Geometry/Shapes.hpp`](./Shapes.hpp)

```cpp
class Point {};
class Line {};
class Segment {};
struct Polygon {};
```

### [`include/Geometry/Mesh.hpp`](./Mesh.hpp)

```cpp
struct Element {};
struct Mesh {};
```

## Methods

### [`include/Geometry/Shapes.hpp`](./Shapes.hpp)

```cpp
// Distances.
Real distance(const Point &, const Point &);
Real distance(const Point &, const Line &);
Real distance(const Point &, const Segment &);

// Intersections.
std::vector<Point> intersections(const Line &, const Line &);
std::vector<Point> intersections(const Line &, const Segment &);
std::vector<Point> intersections(const Line &, const Polygon &);

// Lines.
Line bisector(const Point &, const Point &);
Line normal(const Line &, const Point &);
Line normal(const Segment &, const Point &);

// Polygons.
Polygon collapse(const Polygon &, const Point &);
Polygon collapse(const Polygon &, const Segment &);
std::vector<Point> reflections(const Polygon &, const Point &);
Polygon reduce(const Polygon &, const Line &, const Point &);
```

### [`include/Geometry/Voronoi.hpp`](./Voronoi.hpp)

```cpp
// Voronoi.
std::vector<Polygon> voronoi(const Polygon &, const std::vector<Point> &, const bool &reflect = false);
std::vector<Polygon> voronoi_random(const Polygon &, const std::size_t &, const bool &reflect = false);
std::vector<Polygon> voronoi_uniform(const Polygon &, const std::size_t &, const bool &reflect = false);

// Triangulations.
std::vector<Polygon> triangulate(const Polygon &);
std::vector<Polygon> triangulate(const std::vector<Polygon> &);
```

### [`include/Geometry/Mesh.hpp`](./Mesh.hpp)

```cpp
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
```