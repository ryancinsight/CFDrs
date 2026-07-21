# Geometry Schematics

`cfd-schematics` generates parametric SVG schematics for all key CFDrs
geometry primitives.  Run examples from the `cfd-schematics` crate root.

| Example | Output |
|---------|--------|
| `comprehensive_serpentine_demo` | Full serpentine channel family (aspect, turns, width) |
| `unified_generator_demo` | All topology types from a single generator API |
| `frustum_channel_demo` | Frustum (truncated-cone) Venturi geometry |

```bash
cd crates/cfd-schematics
cargo run --example comprehensive_serpentine_demo
cargo run --example topology/unified_generator_demo
cargo run --example venturi/frustum_channel_demo
```
