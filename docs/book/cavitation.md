# Cavitation: Liquid-Vapour Phase Transition

Cavitation is modeled as a **phase transition** whose rate depends on local
pressure-vs-vapor-pressure:

```rust
pub trait Cavitation {
    fn inception_pressure<F: FloatElement>(&self) -> F;
    fn collapse_rate<F: FloatElement>(&self, p: F) -> F;
}
```

The default model is a **Rayleigh–Plesset** closure; an
**Eulerian–Eulerian** variant is available for strong cloud-cavitation
flows.

Damage accumulation is integrated per timestep as

    damage += (inception_indicator) · (collapse_rate · dt)

and reported through the [`cfd-3d::posterior::damage`] hook.  Outputs are
dimensionally consistent (kg/m²/s) so that downstream fatigue models can
post-process without re-scaling.

## Examples Referenced by This Chapter

- [`venturi_cavitation`](examples/venturi_cavitation.md) — venturi cavitation
  inception and collapse.
- [`simple_cavitation`](examples/simple_cavitation.md) — single-bubble
  Rayleigh–Plesset in a constrained channel.
- [`cavitation_damage_simulation`](examples/cavitation_damage_simulation.md)
  — long-time damage integration under turbulent inflow.
- [`1d_venturi_cavitation`](examples/1d_venturi_cavitation.md) — 1-D
  reduced cavitation model.

## Further Reading

- [`cfd-3d` cavitation module](../../crates/cfd-3d/src/cavitation.rs)
- For the unified `Integrate` trait, turbulence closures, and multiphase
  coupling, see [Turbulence Models and Cavitation](turbulence_multiphase.md).
