# cfd-3d

3D solvers for the [CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

## What it provides

- **`fem`** — finite element discretization, including Stokes flow.
- **`spectral`** — spectral methods built on `apollo-fft` and `apollo-nufft`.
- **`ibm`** — immersed boundary method.
- **`vof`, `level_set`, `multiphase`** — interface tracking: volume of fluid with
  piecewise-linear reconstruction, signed-distance level sets, and phase coupling.
- **`cavitation`** — cavitation modeling and damage estimation.
- **`bifurcation`, `trifurcation`, `serpentine`, `venturi`, `cascade`** — device
  geometries resolved in 3D.
- **`blueprint_integration`** — consumes `cfd-schematic-mesh` output so a schematic
  can be meshed and solved without leaving the suite.

Meshes come from `gaia-mesh` (re-exported through the workspace as `cfd-mesh`).

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
