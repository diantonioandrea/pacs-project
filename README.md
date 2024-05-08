# Advanced Programming for Scientific Computing - Project

_Adaptive HP Discontinuous Galërkin Algorithms_

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [Test Compilation](#test-compilation)

## Introduction

This repository presents an implementation of an adaptive HP Discontinuous Galërkin method.

## Overview

Key components include:

- `include/`: Definitions of methods and classes.
    - `Type.hpp`: Definitions of the custom concepts.
    - `Algebra.hpp`: Methods for matrices and vectors.
    - `Geometry.hpp`: Definitions of geometrical objects and methods.
    - `Matrix.hpp`: Implementation of matrices.
    - `Sparse.hpp`: Implementation of sparse matrices.
    - `Vector.hpp`: Implementation of vectors.
    - `Voronoi.hpp`: Definitions of methods for evaluating the Voronoi diagram 
    - `Mesh.hpp`: Definitions of mesh class and methods.
    - `Quadrature.hpp`: Definitions of quadrature methods and algorithms.
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

- `domain`: Creates a mesh over `[-1, 1] x [-1, 1]` and dumps it in `mesh.poly`, readable by `polyplot.py`.
- `gauss`: Evaluates and prints some **Gauss-Legendre** quadrature nodes and weights over `[-1, 1] x [-1, 1]`.
- `algebra`: Solves a simple linear system using the _Conjugate Gradient Method_.

Tests can be run by:

    ./TEST