# Teal

Teal is a general-purpose computational fluid dynamics library for solving hyperbolic-parabolic
conservation laws on 3D unstructured meshes using the finite-volume method. It includes a simple
mesh generator and integrates with GMSH for complex geometries. Built-in solvers cover the
compressible Euler and Navier-Stokes equations, with support for user-defined source terms. The teal
API is designed for flexibility, making it straightforward to extend to other systems of equations.

## Features

- **Meshes**
    - 3D unstructured mixed meshes
    - Built-in Cartesian mesh generator
    - Readers for GMSH and internal formats
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

## Getting started

Clone and build:

```bash
git clone https://github.com/hheinzer/teal.git
cd teal
make
```

Executables are placed in the `bin` directory.

The best way to get started with teal is by running test cases in [run](run/). For example, the [Sod
shock tube](https://en.wikipedia.org/wiki/Sod_shock_tube):

```bash
mpirun -n 4 bin/riemann/sod
```

The test cases demonstrate all capabilities of teal. If your editor supports "jump to definition,"
you can follow the program flow directly in the source.

Check the [installation](INSTALL.md) instructions if you run into compilation issues.

## Documentation

For now, there is no dedicated documentation. In my opinion, function and variable names together
with inline comments should be enough to understand what's going on. If there is enough interest in
the project, I'll invest time in proper documentation.

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
4. Open a pull request to the main branch

## Acknowledgments

Teal is a successor to [ccfd](https://github.com/hheinzer/ccfd), which itself is a C rewrite of
[cfdfv](https://github.com/flexi-framework/cfdfv). While teal builds on the concepts of these
projects, it does not share any code with them.

The development of teal has been guided by foundational literature in computational fluid dynamics.
Key resources include *Computational Fluid Dynamics: Principles and Applications* by J. Blazek and
*Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction* by E. F. Toro.

## License

Teal is licensed under the GPL-3.0 [license](LICENSE).
