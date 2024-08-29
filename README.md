# Teal

An open source finite volume framework for hyperbolic–parabolic conservation laws.

## Features

- 3D unstructured meshes
- 2nd order accuracy in space and time
- Modular equation system design
- Parallelized with MPI

## Installation

You will require the following dependencies:

- [Make](https://www.gnu.org/software/make/)
- [C compiler](https://clang.llvm.org/)
- [MPI](https://www.open-mpi.org/)
- [GMSH](https://gmsh.info/)
- [METIS](https://github.com/KarypisLab/METIS)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [Python](https://www.python.org/) (optional)
- [ParaView](https://www.paraview.org/) (optional)

Once everything is installed, compile teal with `make` and run simulations with `mpirun`.

## Usage

Since teal is written as a library, simulations are defined in C, not through input files.
Typically, you write a small `main` function to set up and run the simulation. Check the [run](run/)
folder for examples of how to use teal. You can also use command-line arguments or create an input
file parser to parametrize your simulations if needed.

## Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue. To
contribute code:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes and push to your fork.
4. Open a pull request to the main branch.

## Acknowledgments

Teal is a successor to [ccfd](https://github.com/hheinzer/ccfd), which itself is a C rewrite of
[cfdfv](https://github.com/flexi-framework/cfdfv). While teal builds on the concepts from these
projects, it does not share any code with them.

The development of teal has been guided by foundational literature in computational fluid dynamics.
Key resources include *Computational Fluid Dynamics: Principles and Applications* by Jiri Blazek and
*Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction* by E. F. Toro.

## License

Teal is licensed under the GPL-3.0 [license](LICENSE).
