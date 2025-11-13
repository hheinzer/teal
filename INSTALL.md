# Installation

These instructions assume a POSIX environment (Linux or macOS) with MPI. If you want to run teal on
Windows, I recommend that you use WSL.

## Requirements

- C tool chain and build tools: `make` plus either Clang or GCC
- MPI implementation with `mpicc` (Open MPI or MPICH)
- HDF5 **built with MPI support** (headers + libs)
- METIS and ParMETIS (both installed; ParMETIS must match your METIS build)
- libunwind for your platform, adapt the Makefile if it's not x86_64

> For ParMETIS, both `i64=0` and `i64=1` builds are supported. For very large meshes (when 32-bit
> indices overflow), use the `i64=1` build.

## Compilers

The compiler is selected in the Makefile; change it there if needed. The project compiles under both
Clang and GCC, although GCC throws some warnings, which appear bogus to me. It is also possible to
use TCC, albeit with a performance hit.
