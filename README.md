# Advanced Programming for Scientific Computing - Project

_Adaptive HP Discontinuous Galërkin Algorithms_

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [Test Compilation](#test-compilation)
- [Usage](#usage)
    - [Using **polyplot**](#using-polyplot)

## Introduction

This repository presents an implementation of an adaptive HP Discontinuous Galërkin method.

## Overview

Key components include:

- `include/`: Definitions of methods and classes.
    - `Type.hpp`: Definitions of the custom concepts and types.
    - `Geometry.hpp`: Definitions of geometrical objects and methods.
    - `Matrix.hpp`: Implementation of matrices.
    - `Sparse.hpp`: Implementation of sparse matrices.
    - `Vector.hpp`: Implementation of vectors.
    - `Voronoi.hpp`: Definitions of methods for evaluating the Voronoi diagram 
    - `Mesh.hpp`: Definitions of mesh class and methods.
    - `Quadrature.hpp`: Definitions of quadrature methods and algorithms.
    - `Basis.hpp`: Definitions of basis functions methods.
    - `Legendre.hpp`: Definitions of Legendre polynomials methods.
    - `Laplacian.hpp`: Definitions of methods for the laplacian matrix.
    - `Forcing.hpp`: Definitions of methods for the forcing term.
    - `Penalty.hpp`: Definitions of methods for the penalty coefficients.
    - `Functor.hpp`: Definitions of the functor class.
    - `Solution.hpp`: Definitions of methods for working with the numerical solution.
    - `Errors.hpp`: Definitions of methods for evaluating the errors.
- `src/`: Implementations of the methods and classes defined under `include/`.
- `scripts/`: Useful scripts.
    - `polyplot.py`: Plots a mesh dumped with `Mesh::write(const std::string &filename)`
- `test/`: Some tests for the code. See [Test Compilation](#test-compilation).

## Setup

### Cloning the Repository

To begin, clone the repository from [here](https://github.com/diantonioandrea/pacs-project):

    git clone git@github.com:diantonioandrea/pacs-project.git

### Test Compilation

All tests can be compiled with:

    make test

Some meaningful tests may include:

- `test_domain`: Creates a mesh over `[-1, 1] x [-1, 1]` and dumps it in `mesh.poly`.
- `test_gauss`: Evaluates and prints some **Gauss-Legendre** quadrature nodes and weights over `[-1, 1] x [-1, 1]`.
- `test_algebra`: Solves a simple linear system `Ax = b`.

Single tests can be run by:

    ./executables/test_TESTNAME

while the whole set of tests can be run by:

    make run

## Usage

### Using **polyplot**

Mesh files can be plotted using **polyplot** by:

    ./scripts/polyplot.py file.poly