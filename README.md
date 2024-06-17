# Advanced Programming for Scientific Computing - Project

_Adaptive HP Discontinuous Galërkin Algorithms_

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)

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