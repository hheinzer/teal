# Installation

These instructions assume a POSIX environment (Linux or macOS) with MPI. If you want to run teal on
Windows, I recommend that you use WSL.

## Requirements

- C tool chain and build tools: `make` plus either Clang or GCC
- MPI implementation with `mpicc` (Open MPI or MPICH)
- HDF5 **built with MPI support** (headers + libs)
- METIS and ParMETIS (both installed; ParMETIS must match your METIS build)
- Gmsh (to generate meshes from `.geo` in `run/`)
- Python + matplotlib/pyvista/mpi4py (for `tools/plot` and helper scripts)

## Compilers

The compiler is selected in the Makefile; change it there if needed. The project compiles under both
Clang and GCC, although GCC throws some warnings, which appear bogus to me. It is also possible to
use TCC, albeit with a performance hit.
