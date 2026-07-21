# 3-D Flows

Three-dimensional CFD examples live in `cfd-3d`.  They exercise the full
3-D pressure-velocity coupling, turbulence, and cavitation pathways.

| Example | Physics |
|---------|---------|
| `spectral_poisson_3d` | Spectral Poisson solver on a 3D periodic domain |
| `bifurcation_3d_blood` | Non-Newtonian blood in a symmetric 3-D bifurcation |
| `venturi_3d_cavitation` | Cavitation inception in a 3-D Venturi constriction |
| `serpentine_3d_dean` | Dean-vortex mixing in a 3-D serpentine channel |

```bash
cargo run -p cfd-3d --example spectral_poisson_3d
cargo run -p cfd-3d --example bifurcation_3d_blood
cargo run -p cfd-3d --example venturi_3d_cavitation
cargo run -p cfd-3d --example serpentine_3d_dean
```
