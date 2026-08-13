# cfd-schematics

Schematic network design for microfluidic and millifluidic devices, part of the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

A `NetworkBlueprint` is the single source of truth for device topology and geometry.
`cfd-1d`, `cfd-2d`, and `cfd-schematic-mesh` all consume it, so a device is described
once and then solved at whichever fidelity the question needs.

## What it provides

- **`topology`** — nodes, channels, and network specification types.
- **`domain`, `application`, `infrastructure`, `interface`** — the layered blueprint
  model, its builders, and its persistence boundary.
- **`geometry`** — series, parallel, and selective-tree geometry generators.
- **`config`** — configuration manifests and physical constants.
- **`heatmap`, `visualizations`** — field rendering over the schematic. Normalized
  color laws come from `iris`; this crate keeps the CFD field semantics and sparse
  node/channel association, borrowing solver maps through `Cow` and precomputing
  finite scalar ranges before rendering.
- **`state_management`** — schematic editing state.

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
