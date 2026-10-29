# advdiff2d-hpc-lab

[![CI](https://github.com/ldenarie/advdiff2d-hpc-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/ldenarie/advdiff2d-hpc-lab/actions/workflows/ci.yml)

> **Work in progress — learning project.**

A 2D advection-diffusion kernel in C++, used to learn and practise parallel
programming (OpenMP, then MPI), vectorization and performance measurement.

## What this is

- A small, tested explicit finite-difference solver for a passive tracer on a
  periodic 2D grid: first-order upwind advection, centred diffusion, forward
  Euler. See [docs/model.md](docs/model.md).
- A playground for parallelization, with reproducible measurements and an
  honest discussion of their limits.

## What this is not

- Not a weather or climate model, nor a simplified version of one.
- Not a large-scale scalability study: measurements run on small free cloud
  machines (a few cores, one node).

## Status

| Step | Status |
|---|---|
| Sequential kernel and tests | in progress |
| Profiling | planned |
| OpenMP | planned |
| Vectorization | planned |
| MPI (domain decomposition) | planned |

## Build and test

Requires CMake >= 3.22, a C++17 compiler and Catch2 v3 (`catch2` on
Ubuntu). The repository ships a dev container with everything installed.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/advdiff2d nx=512 ny=512 steps=200
```

## Use of AI assistance

The project skeleton and tooling (CMake, CI, dev container, specification
draft) were prepared with the help of an AI assistant. The kernels,
parallelization, tests and analysis are written by the author.

## License

MIT, see [LICENSE](LICENSE).
