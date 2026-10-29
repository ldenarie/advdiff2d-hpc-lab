# Model and numerical scheme

This document is the specification of the exercise. The code must follow it;
any deviation is documented here first.

## 1. Equation

Transport of a passive scalar (tracer) `c(x, y, t)` by a uniform, constant
velocity `(u, v)` with isotropic diffusion `kappa`:

```
dc/dt + u dc/dx + v dc/dy = kappa (d2c/dx2 + d2c/dy2)
```

- Domain: `[0, Lx) x [0, Ly)`, **periodic** in both directions.
- `u`, `v` may be positive, negative or zero. `kappa >= 0`.
- The velocity is uniform on purpose: it keeps the kernel a pure stencil
  (memory-bound), which is the point of the exercise. A variable velocity field
  (e.g. solid-body rotation) is out of scope for v0.x.

## 2. Discretization

### Grid

- `nx x ny` cells, `dx = Lx / nx`, `dy = Ly / ny`.
- Values are cell-centred: `x_i = (i + 0.5) dx`, `y_j = (j + 0.5) dy`,
  with `i in [0, nx)`, `j in [0, ny)`.

### Time stepping

Forward Euler, explicit:

```
c_new(i,j) = c(i,j) + dt * ( A(i,j) + D(i,j) )
```

### Advection: first-order upwind

```
A_x(i,j) = -u * (c(i,j)   - c(i-1,j)) / dx   if u >= 0
A_x(i,j) = -u * (c(i+1,j) - c(i,j))   / dx   if u <  0
```

and the same in `y` with `v`, `dy`, `j`. `A = A_x + A_y`.

### Diffusion: centred 5-point Laplacian

```
D(i,j) = kappa * ( (c(i+1,j) - 2 c(i,j) + c(i-1,j)) / dx^2
                 + (c(i,j+1) - 2 c(i,j) + c(i,j-1)) / dy^2 )
```

## 3. Stability and maximum principle

Define the dimensionless numbers

```
Cx = |u| dt / dx        Cy = |v| dt / dy
Dx = kappa dt / dx^2    Dy = kappa dt / dy^2
S  = Cx + Cy + 2 Dx + 2 Dy          (stability number)
```

Rewriting the update shows that `c_new(i,j)` is a weighted sum of `c(i,j)` and
its four neighbours, with weights that sum to 1. The weight of `c(i,j)` is
`1 - S`; all the others are non-negative. Therefore:

- if `S <= 1`, the update is a **convex combination**: the scheme is stable,
  and the new values stay within `[min c, max c]` (discrete maximum principle);
- the weights sum to 1 and the domain is periodic, so the **total mass**
  `sum c(i,j) dx dy` is conserved exactly, up to floating-point rounding.

The driver chooses `dt` from a target `S` (default 0.9).

## 4. Data layout

- One layer of **ghost (halo) cells** around the interior: storage is
  `(nx + 2) x (ny + 2)` doubles.
- Row-major, `x` is the fastest index: interior cell `(i, j)` is at
  `(j + 1) * (nx + 2) + (i + 1)`.
- Before each step, `apply_periodic_bc` copies the opposite interior edges into
  the ghost cells. Corners are not needed (5-point stencil).
- Two buffers (`c`, `c_new`), swapped after each step (double buffering).

Why ghost cells rather than `% nx` in the inner loop: the inner loop stays
branch-free and contiguous (vectorizable), and the MPI version will replace the
periodic copy in `y` by a halo exchange between neighbouring ranks without
touching the kernel.

## 5. Verification (what the tests must check)

| Test | Setup | Expected |
|---|---|---|
| Uniform field | constant `c`, any `u, v, kappa`, `S <= 1` | unchanged, exactly |
| Periodic BC | known interior values | each ghost cell equals the opposite interior cell |
| Mass conservation | Gaussian, `u, v, kappa` all non-zero, `S = 1`, 500 steps | relative mass error `< 1e-12` |
| Maximum principle | same run | `min` and `max` stay within the initial `[min, max]` (tolerance `1e-14`) |
| Diffusion order | `u = v = 0`, Fourier mode `1 + sin(2 pi x) sin(2 pi y)`, `dt = 0.2 dx^2 / kappa`, grids 32, 64, 128 | observed order of the RMS error vs. exact decay `exp(-8 pi^2 kappa t)`: between 1.8 and 2.2 |
| Advection order | `kappa = 0`, same mode, `dt` from `S = 0.5`, grids 64, 128, 256 | observed order vs. exact translated mode: between 0.9 and 1.1 |

Choose `dt` so that the final time is reached in an integer number of steps
(compute the number of steps, then adjust `dt = T / steps`).

The observed order between two grids `n` and `2n` is
`log2(error_n / error_2n)`.

## 6. Out of scope

This is not a weather model: no variable wind, no physics, no
semi-Lagrangian or spectral method, no vertical dimension. The kernel is a
representative memory-bound stencil, chosen to practise parallelization and
performance measurement.
