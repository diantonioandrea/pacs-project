# Advanced Programming for Scientific Computing - Project

_The hp-Adaptive Discontinuous Galërkin Method_

## Introduction

This repository provides an implementation of the _hp-adaptive_ discontinuous Galërkin method for the Poisson problem, serving as the foundation for the _PacsHPDG_ library.

:warning: Refer to the [report](#compiling-the-report) for a detailed explanation of the background and results related to this project.

## Table of Contents

- [Introduction](#introduction)
- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [**Compilation and Execution**](#compilation-and-execution)
        - [Compilation Flags](#compilation-flags)
    - [Compiling the Code into a Library](#compiling-the-code-into-a-library)
    - [**Compiling the Report**](#compiling-the-report)
    - [Compiling the Presentation](#compiling-the-presentation)
- [Using the Reporitory](#using-the-repository)
    - [Domains](#domains)
    - [**Examples**](#examples)
    - [Scripts](#scripts)
- [Using the Code](#using-the-code)
    - [Generating a Mesh](#generating-a-mesh)
    - [Solving the Poisson Problem](#solving-the-poisson-problem)
    - [Mesh Refinement](#mesh-refinement)
- [**Notes to the Reader**](#notes-to-the-reader)
    - [On the Implementation of Basic Objects](#on-the-implementation-of-basic-objects)
    - [On the Adaptation from **lymph**](#on-the-adaptation-from-lymph)
    - [On the Examples Structure](#on-the-examples-structure)
    - [On the Custom Laplacian Solver](#on-the-custom-laplacian-solver)
    - [On Parallelism](#on-parallelism)

:warning: Make sure to take a look at [Notes to the Reader](#notes-to-the-reader) as they provide insight into some design choices about the code.

## Overview

The key components are as follows:

- `include/`: Holds definitions for the structures and methods utilized in the repository.
    - [`include/PacsHPDG/Algebra/`](./include/PacsHPDG/Algebra/): Structures and methods for vectors, matrices and linear solvers.
    - [`include/PacsHPDG/Geometry/`](./include/PacsHPDG/Geometry/): Tools for working with polygons and meshes.
    - [`include/PacsHPDG/Fem/`](./include/PacsHPDG/Fem/): Finite element structures and methods.
    - [`include/PacsHPDG/Laplacian/`](./include/PacsHPDG/Laplacian/): Implementation details for the Poisson problem.
    - [`include/PacsHPDG/Statistics/`](./include/PacsHPDG/Statistics/): Statistics related tools.
- `src/`: Contains the primary implementations for the repository’s structures and methods.
- `data/`: Includes sample meshes for simple domains.
- `domains/`: Stores scripts for generating sample meshes.
- `examples/`: Provides the main examples for using the repository.
- `snippets/`: Simpler examples for the report.
- `test/`: Contains scripts for testing fundamental features.
- `scripts/`: Contains Python scripts for meshes, solutions and errors visualization.
- `templates/`: Contains TikZ templates.
- `report/`: Contains the LaTeX report on the project.
- `presentation/`: Contains the LaTeX presentation on the project.

Every directory under `include/PacsHPDG/` has a `README.md` that lists the classes, structures, and methods introduced in that particular category.

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/pacs-project):

```bash
git clone git@github.com:diantonioandrea/pacs-project.git
```

### Compilation and Execution

The code is written to work with the C++ standard library. It has an optional dependency on [_OpenMP_](https://www.openmp.org). The `Makefile` is designed to check for the presence of GCC, modules or a custom installation of _OpenMP_ using the definition of `$(OpenMP)` set to `/path/to/libomp`[^clang].

[^clang]: Tested on Apple's clang.

Compile everything[^compilation] with:

[^compilation]: The `-j` option is advised.

```bash
make
```

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

#### Compilation Flags

:warning: Be sure to review the [Makefile](./Makefile) for additional options to customize the code's behavior.

The code uses the following custom compilation flags:

- `-DNVERBOSE`: Disables verbose output.
- `-DNSOLUTIONS`: Disables solution output.

### Compiling the Code into a Library

The code can also be compiled into the static library `PacsHPDG` under `./lib` by:

```bash
make lib
```

which can be installed locally under `~` by:

```bash
make install
```

### Compiling the Report

The report for the project can be compiled by:

```bash
cd report
make
```

Ensure that [_latexmk_](https://ctan.org/pkg/latexmk/) is installed.

### Compiling the Presentation

Similarly, the presentation for the project can be compiled by:

```bash
cd presentation
make
```

Again, ensure that _latexmk_ is installed.

## Using the Repository

Here is a brief guide on using the executables in this repository. Note that tests were written solely to verify specific parts of the code, so their usage is not covered here.

All output generated by the following codes is located in the `output/` directory.

### Domains

Codes written under `domains/` have the purpose of generating meshes over the domains used for testing the code, specifically the square and L-shaped domains.

This repository provides the following:

- `square_domain.cpp`
- `lshape_domain.cpp`

For example, usage would look like:

```bash
./executables/square_domain.out 250
```

This command generates a square mesh over $[0, 1] \times [0, 1]$ with $N = 250$.

```bash
./executables/lshape_domain.out 125
```

This command generates an L-shaped mesh over $[-1, 1] \times [-1, 1] \setminus [0, 1] \times [-1, 0]$ with $N = 125$.

### Examples

Examples are divided into the following categories based on their domains and meshes:

1. Uniform meshes:
    - `square_smooth.cpp` Square domain, smooth solution[^solutions].
    - `square.cpp` Square domain, non-smooth solution.
    - `lshape_smooth.cpp` L-shaped domain, smooth solution.
    - `lshape.cpp` L-shaped domain, non-smooth solution.
2. _h-adaptively_ refined meshes with _a priori_ estimates based on $L^2$ error:
    - `square_h.cpp` Square domain, non-smooth solution.
    - `lshape_h.cpp` L-shaped domain, non-smooth solution.
3. _h-adaptively_ refined meshes with _a priori_ estimates based on $H^1$ error:
    - `square_gh.cpp` Square domain, non-smooth solution.
    - `lshape_gh.cpp` L-shaped domain, non-smooth solution.
4. _h-adaptively_ refined meshes with _a posteriori_ estimates:
    - `square_eh.cpp` Square domain, non-smooth solution.
    - `lshape_eh.cpp` L-shaped domain, non-smooth solution.
5. _hp-adaptively_ refined meshes with _a posteriori_ estimates:[^hp]
    - `square_hp.cpp` Square domain, non-smooth solution.
    - `lshape_hp.cpp` L-shaped domain, non-smooth solution.

[^solutions]: Smooth and non-smooth solutions are discussed in the report.

[^hp]: The polynomial degree specified for _hp-adaptively_ refined meshes is treated as the starting degree.

Category _1_ requires the user to specify the polynomial degree. An example command is:

```bash
./executables/square.out 3
```

which solves the Poisson problem on a sequence of uniformly refined square meshes with $k = 3$.

Categories _2, ..., 5_ require the user to specify the polynomial degree and optionally a starting mesh identified by its elements. Meshes are stored under `data/square/` or `data/lshape/`. An example command is:

```bash
./executables/square_h.out 3 250
```

which solves the Poisson problem on a sequence of _h-adaptively_ refined square meshes with $k = 3$ and $N_0 = 250$.

:warning: These examples contribute to the graphs presented in the report.

### Scripts

This repository includes some Python scripts to help visualize error trends, meshes and solutions.

Scripts are divided into the following categories based on their function:

1. Mesh-related:
    - `polyplot.py`: Requires a `.poly` file from which it plots a mesh. Accepts `--degrees` for _hp-adaptively_ refined meshes.
2. Solution-related:
    - `solplot.py`: Requires a `.sol` file from which it plots a solution.
3. Error trends:
    - `errorvsize.py`: Requires *one* `.error` file from which it plots error trends versus mesh size.
    - `errorvdofs.py`: Requires *one or two* `.error` files from which it plots error trends versus DOFs.
    - `hpvdofs.py`: Requires *one or two* `.error` files from which it plots error trends versus DOFs. Works for _hp-adaptively_ refined meshes.
4. TikZ wrappers:
    - `evs_tikz.py`: `errorvsize.py` TikZ wrapper.
    - `evd_tikz.py`: `errorvdofs.py` TikZ wrapper.
    - `hpvd_tikz.py`: `hpvdofs.py` TikZ wrapper.

## Using the Code

To use the code, include `<PacsHPDG.hpp>` which provides all necessary tools.

### Mesh Generation

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

The diagram can be retrieved using the same `mesh_diagram` function:

```cpp
std::vector<pacs::Polygon> retrieved = pacs::mesh_diagram("output/square.poly");
```

For non-convex domains, enable point reflection:

```cpp
std::vector<pacs::Polygon> diagram = pacs::mesh_diagram(domain, 100, reflect=true);
```

### Solving the Poisson Problem

First, build the Laplacian matrix from your mesh[^laplacian][^real]:

[^laplacian]: `laplacian` also builds some other matrices used for error evaluation.
[^real]: `pacs::Real` wraps `long double`.

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

### Mesh refinement

After solving the Poisson problem, you can evaluate _a posteriori_ error estimates, allowing for adaptive mesh refinement.

First, evaluate the estimates:

```cpp
pacs::Estimator estimator{mesh, M, numerical, source, dirichlet, {dirichlet_x, dirichlet_y}};
```

Pass the current mesh, the mass matrix, the source, and the Dirichlet boundary condition along with its derivatives to the `pacs::Estimator` constructor.

Use the `mesh_refine` function to _hp-refine_ the mesh:

```cpp
mesh_refine(mesh, estimator);
```

## Notes to the Reader

### On the Implementation of Basic Objects

This repository implements basic objects such as vectors and matrices, which are also found in many other libraries. Additionally, it includes a polygonal mesher, though it may not be as powerful as existing alternatives. The decision to develop these components from scratch was driven by a desire to minimize dependencies and to create a project that I could truly call my own, showcasing what I learned during the course.

### On the Adaptation from **lymph**

The implementation of the Laplacian problem and FEM tools in this project was adapted from the [**lymph**](https://lymph.bitbucket.io) library.

### On the Examples Structure

All the examples have the same structure, and they may seem quite repetitive and superfluous. There is no a priori reason to treat every single example as a different file, given the source code similarities despite the different functions. However, this approach makes each example easy to comprehend and does not pose an issue since examples do not need to be maintained, which again wouldn't be difficult thanks to modern IDEs and editors' capabilities.

### On the Custom Laplacian Solver

All the adaptive examples use `pacs::lapsolver` instead of a manual `pacs::solve`. This is a wrapper for `GMRES` with a `DBI`[^dbi] preconditioner, which is particularly useful given the ill-conditioning that _hp-adaptive_ methods suffer from.

[^dbi]: A generic name for a Block-Jacobi preconditioner.

### On Parallelism

Parallelism through the STL plays a secondary role, primarily utilized for common operations involving standard containers. In contrast, parallelism via _OpenMP_ is fundamental, significantly boosting the performance of the polygonal mesher, the `DB` solver, and the `DBI` preconditioner.