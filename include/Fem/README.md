# `include/Laplacian/`

## Introduction

This directory contains the implementation of FEM-related components for the project.

## Classes and structs

### [`include/Fem/Functor.hpp`](./Functor.hpp)

```cpp
class Functor {};
class TwoFunctor {};
```

### [`include/Fem/Solution.hpp`](./Solution.hpp)

```cpp
struct Solution {};
```

## Methods

### [`include/Fem/Functor.hpp`](./Functor.hpp)

```cpp
inline Real zero(const Real &, const Real &) { return 0.0; }
```

### [`include/Fem/Basis.hpp`](./Basis.hpp)

```cpp
std::array<Matrix<Real>, 3> basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
Matrix<Real> lap_basis_2d(const Mesh &, const std::size_t &, const std::array<Vector<Real>, 2> &);
```

### [`include/Fem/Legendre.hpp`](./Legendre.hpp)

```cpp
Vector<Real> legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> grad_legendre(const Vector<Real> &, const std::size_t &);
Vector<Real> lap_legendre(const Vector<Real> &, const std::size_t &);
```

### [`include/Fem/Quadrature.hpp`](./Quadrature.hpp)

```cpp
// Nodes.
std::array<Vector<Real>, 2> gauss_legendre(const Real &, const Real &, const std::size_t &);

// Reference nodes.
std::array<Vector<Real>, 2> quadrature_1d(const std::size_t &);
std::array<Vector<Real>, 3> quadrature_2d(const std::size_t &);
```

### [`include/Fem/Modal.hpp`](./Modal.hpp)

```cpp
Vector<Real> modal(const Mesh &, const Functor &);
```

### [`include/Fem/Penalty.hpp`](./Penalty.hpp)

```cpp
Vector<Real> penalty(const Mesh &, const std::size_t &);
```