# Advanced Programming for Scientific Computing - Project

_Adaptive HP Discontinuous Galërkin Algorithms_

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
- [Usage](#usage)
 - [Generating a Mesh](#generating-a-mesh)

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

## Setup

### Cloning the Repository

To begin, clone the repository from [here](https://github.com/diantonioandrea/pacs-project):

    git clone git@github.com:diantonioandrea/pacs-project.git

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

Save the generated mesh to a `.poly` file:

```cpp
pacs::Mesh mesh{domain, diagram};
mesh.write("output/square.poly");
```

For non-convex domains, enable point reflection:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100, reflect=true);
```