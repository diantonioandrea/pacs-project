# `include/Laplacian/`

## Introduction

This directory contains the implementation of algebra-related components for the project.

## Classes and structs

### [`include/Laplacian/Errors.hpp`](./Errors.hpp)

```cpp
struct Error {};
```

### [`include/Laplacian/Estimators.hpp`](./Estimators.hpp)

```cpp
struct Estimator {};
```

## Methods

### [`include/Laplacian/Laplacian.hpp`](./Laplacian.hpp)

```cpp
// Matrices.
std::array<Sparse<Real>, 3> laplacian(const Mesh &);

// Blocks.
std::vector<std::array<std::vector<std::size_t>, 2>> block_mass(const Mesh &);
```

### [`include/Laplacian/Forcing.hpp`](./Forcing.hpp)

```cpp
Vector<Real> forcing(const Mesh &, const Functor &, const Functor &dirichlet = Functor{});
```

### [`include/Laplacian/Refine.hpp`](./Refine.hpp)

```cpp
void mesh_refine(Mesh &, const Estimator &, const Real &refine = 0.75, const Real &speed = 1.0);
```

### [`include/Laplacian/Solvers.hpp`](./Solvers.hpp)

```cpp
Vector<Real> lapsolver(const Mesh &, const Sparse<Real> &, const Vector<Real> &, const Real &TOL = 1E-15);
```