![](.teal.svg)

Teal is a general-purpose computational fluid dynamics library for solving hyperbolic-parabolic
conservation laws on 3D unstructured meshes using the finite-volume method. It ships with a simple
mesh generator, integrates with Gmsh for complex geometries, and includes solvers for compressible
Euler and Navier-Stokes equations with support for user-defined source terms. The teal API is
designed for flexibility, making it straightforward to extend to other systems of equations.

## Highlights

- **Meshes**
    - 3D unstructured mixed meshes
    - Built-in Cartesian mesh generator
    - Parallel reader for Gmsh
    - Extensible to additional formats
- **Accuracy**
    - 2nd-order spatial accuracy
        - Least-squares gradient reconstruction
    - Up to 3rd-order temporal accuracy
        - Low-storage explicit Runge-Kutta schemes
        - Implicit Euler with Newton-GMRES solver
- **Equations**
    - Compressible Euler and Navier-Stokes
    - Modular design for adding other conservation laws

## Performance

Teal is a toy code that I write for *fun*. Nevertheless, this question is bound to come up
eventually. To test the [strong scaling](https://hpc-wiki.info/hpc/Scaling) of teal, I run the
[Taylor-Green vortex](run/taylor_green_vortex/) test case with different resolutions and vary the
number of MPI ranks, keeping the problem size fixed in each case. For each resolution, speedup and
efficiency are computed relative to the 128-rank run.

| cells | ranks | time      | speedup | efficiency [%] |
|-------|-------|-----------|---------|----------------|
| 128^3 | 128   | 1m 18.0s  | 1.00    | 100.0          |
| 128^3 | 256   | 39.8s     | 1.96    | 98.0           |
| 128^3 | 512   | 23.9s     | 3.26    | 81.4           |
| 128^3 | 1024  | 26.7s     | 2.92    | 36.6           |
| 256^3 | 128   | 20m 19.4s | 1.00    | 100.0          |
| 256^3 | 256   | 10m 44.2s | 1.89    | 94.7           |
| 256^3 | 512   | 5m 33.7s  | 3.65    | 91.3           |
| 256^3 | 1024  | 3m 5.7s   | 6.57    | 82.1           |

## Requirements

- C tool chain and build tools: `make` plus either Clang or GCC (set `CC` in `Makefile`)
- MPI implementation with `mpicc` (Open MPI or MPICH) in your `PATH`
- HDF5 **built with MPI support** (headers + libs)
- METIS and ParMETIS (both installed; ParMETIS must match your METIS build)
- Gmsh (to generate meshes from `.geo` files)

## Getting started

Clone and build:

```bash
git clone https://github.com/hheinzer/teal.git
cd teal
make
```

Executables are placed in the `bin` directory.

The best way to get started with teal is by running test cases in [run](run/). For example:

```bash
mpirun -n 4 bin/supersonic_wedge/run
```

Visualize the resulting VTKHDF files directly in ParaView.

The test cases demonstrate all capabilities of teal. If your editor supports "jump to definition",
you can follow the program flow directly in the source.

## Documentation

For now, there is no dedicated documentation. In my opinion, function and variable names together
with inline comments should be enough to understand what's going on. If there is enough interest in
the project, I'll invest time into proper documentation.

If you want to explore the code, the [src](src/) directory contains the main modules of teal. Each
module's header file provides a brief overview and lists public functions. The folder named after
the module contains the implementation. One exception is the `teal` module which contains the core
tools used across all modules.

## Contributing

Contributions are welcome!

If you find a bug or have a feature request, please open an issue. To contribute code:

1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Commit your changes and push to your fork
4. Open a pull request to the `next` branch

## Acknowledgments

Teal is a successor to [ccfd](https://github.com/hheinzer/ccfd), which itself is a C rewrite of
[cfdfv](https://github.com/flexi-framework/cfdfv). While teal builds on the concepts of these
projects, it does not share any code with them.

The development of teal has been guided by foundational literature in computational fluid dynamics.
Key resources include *Computational Fluid Dynamics: Principles and Applications* by J. Blazek and
*Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction* by E. F. Toro.

## License

Teal is licensed under the GPL-3.0 [license](LICENSE).
