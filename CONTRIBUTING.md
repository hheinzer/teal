# Contributing

Thanks for your interest in teal! Contributions of all kinds help make the solver better, and I’m
glad you’re here to improve it.

## Goals

- Keep teal simple, deterministic, and reproducible across MPI ranks
- Favor clarity over cleverness; descriptive names >> micro-optimizations
- Treat example cases in `run/` as integration tests and documentation
- Prefer changes that improve robustness, add coverage, or clarify APIs

## Environment

- See [INSTALL](INSTALL.md) for required dependencies
- Default build uses Clang via `mpicc` with debug flags; switch to release flags in the
  [Makefile](Makefile) when benchmarking
- Build everything with `make` at repo root

## Workflow

- Develop on a feature branch; keep commits small and focused
- Prefer adding or updating a driver in `run/` to demonstrate new behavior
- Keep PRs and commits scoped; describe what changed and why in the message

## Coding style

- C99, headers in `src/<module>.h`, implementations in `src/<module>/*.c`
- Keep interfaces minimal; pass raw pointers/arrays and keep allocations explicit
- Use existing helpers: arena allocator for transient buffers, `println/verbose/error`, `array_*`
  for reductions, `sync_*` for MPI collectives, ...
- Follow existing naming patterns: `equations_*`, `mesh_*`, `simulation_*`, `euler_*`, ...
- Stick to these types: `long` for indices and sizes, `scalar` for floating point numbers, `vector`
  for xyz scalars, `tuple` for xyz longs; use other types only when matching an external API
- Add brief comments only where intent is non-obvious; avoid redundant narration
- Use `make format` (clang-format) and `make tidy` (clang-tidy) as needed
- Log with `println/verbose`, abort with `error`, avoid stray `printf`
- Align new file/directory names with existing modules; keep public declarations in headers

## Testing and validation

- Run `make` to compile; fix all warnings
- For CFD changes, run relevant test cases before and after the change and mention the results
- When touching I/O, verify VTKHDF outputs still open in ParaView
- Add new regression drivers under `run/` when feasible and keep them small enough for quick runs

## Boundary and physics hooks

- Boundary conditions, convective/viscous fluxes, limiters, and source terms are selected by string
  name; register new ones through the selector plumbing in `equations_*`/`euler_*`
- Keep new hooks MPI-safe and deterministic; avoid hidden global state
- Respect conserved-before-primitive ordering when adding variables or flux terms

## Robots

- Read headers first; they define the contract and naming
- Use `rg` for search and `apply_patch` for scoped edits; do not rewrite unrelated code
- Preserve sanitizer flags and warning levels; avoid adding non-ASCII
- Keep bullet lines terminal-punctuation-free as in this file
- Validate commands are runnable in the current sandbox; prefer fast, local checks
- Mirror the existing logging and allocation patterns; avoid inventing new helpers without need
