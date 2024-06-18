# Advanced Programming for Scientific Computing - Project

_Adaptive HP Discontinuous Galërkin Algorithms_

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [Compilation and Execution](#compilation-and-execution)
- [Usage](#usage)
 - [Generating a Mesh](#generating-a-mesh)
 - [Solving the Poisson Problem](#solving-the-poisson-problem)

## Introduction

This repository presents an implementation of an adaptive HP Discontinuous Galërkin method.

## Overview

The key directories are as follows:

- `include/`: Holds definitions for the structures and methods utilized in the repository.
    - `include/Algebra/`: Structures and methods for vectors, matrices and linear solvers.
    - `include/Geometry/`: Tools for working with polygons and meshes.
    - `include/Fem/`: Finite element structures and methods.
    - `include/Laplacian/`: Implementation details for the Poisson problem.
- `src/`: Contains the primary implementations for the repository’s structures and methods.
- `data/`: Includes sample meshes for simple domains.
- `domains/`: Stores scripts for generating sample meshes.
- `example/`: Provides the main examples for using the repository.
- `test/`: Contains scripts for testing fundamental features.
- `scripts/`: Contains Python scripts for meshes, solutions and errors visualization.

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/pacs-project):

```bash
git clone git@github.com:diantonioandrea/pacs-project.git
```

### Compilation and Execution

Compile examples with:

```bash
make examples
```

Compile tests with:

```bash
make tests
```

Compile mesh generation scripts with:

```bash
make domains
```

Executables are located in `executables/` and their outputs in `output/`.

## Usage

### Mesh Generation

To generate a mesh, include `<Geometry.hpp>` which provides all necessary tools.

Start by defining a domain using points:

```cpp
pacs::Point a{0.0, 0.0};
pacs::Point b{1.0, 0.0};
pacs::Point c{1.0, 1.0};
pacs::Point d{0.0, 1.0};

pacs::Polygon domain{{a, b, c, d}};
```

Next, create a diagram with 100 elements over the square domain:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100);
```

Generate a mesh and save it to a `.poly` file:

```cpp
pacs::Mesh mesh{domain, diagram};
mesh.write("output/square.poly");
```

For non-convex domains, enable point reflection:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100, reflect=true);
```

### Solving the Poisson Problem

To solve the Poisson problem, ensure you have included `<Fem.hpp>` and `<Laplacian.hpp>` for necessary functionalities.

First, build the Laplacian matrix from your mesh:

```cpp
std::array<pacs::Matrix<pacs::Real>, 3> matrices = pacs::laplacian(mesh);
pacs::Matrix<pacs::Real> laplacian = matrices[1];
```

Next, construct the forcing term using specified source and Dirichlet boundary conditions:

```cpp
pacs::Vector<pacs::Real> forcing = pacs::forcing(mesh, source, dirichlet);
```

Here, `source` and `dirichlet` are functions representing the source term and Dirichlet boundary conditions, respectively, which follow the following scheme:

```cpp
pacs::Real function(const pacs::Real &, const pacs::Real &);
```

Finally, solve the linear system to find the solution vector:

```cpp
pacs::Vector<pacs::Real> solution = pacs::solve(laplacian, forcing);
```

This `solution` vector now contains the computed solution to the Poisson problem on the given mesh with specified boundary conditions and source term.