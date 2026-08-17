# cfd-schematic-mesh

The schematic-to-mesh bridge in the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite: it converts a
`cfd-schematics::NetworkBlueprint` into a watertight surface mesh that `cfd-3d` can
solve on.

## What it provides

- **`blueprint_mesh`** — blueprint to watertight surface mesh conversion.
- **`shell_mesh`** — shell/wall generation around the fluid volume.
- **`topology`** — mesh topology construction and validation.
- **`region_map`** — maps mesh regions back to their originating blueprint entities,
  so solver results can be attributed to named channels and nodes.
- **`constraint`** — geometric constraints applied during meshing.
- **`well_plate`** — standard well-plate layouts.
- **`scheme_io`** — schematic import/export.

Meshes are produced through `gaia-mesh` (the workspace `cfd-mesh` dependency).

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
