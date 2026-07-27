# Cell-Separation Force Metrics

## Contract

The public cell-separation family carries dimensional values as Aequitas
quantities:

- `lateral_position` is `Length`;
- `residual_force` is `Force`;
- `dean_drag` is `Force`.

Direct margination and cell-interaction functions use `MassDensity`,
`DynamicViscosity`, `Velocity`, and `Length` for physical inputs. The
three-population entry points use the same contracts plus
`VolumetricFlowRate`; normalized positions, fractions, Reynolds number, and
Dean number remain dimensionless scalars.

- Fahraeus, Fahraeus–Lindqvist, cell-free-layer, and Quemada APIs use typed
  diameter, viscosity, and reciprocal-time inputs and return typed widths or
  dynamic viscosities.
- Plasma-skimming and cross-junction cascade APIs use typed diameters,
  channel geometry, and volumetric flow rates.

`x_tilde_eq`, Reynolds number, Dean number, and separation fractions remain
dimensionless model values. `CellSeparationModel` stores channel width, height,
length, and optional bend radius as `Length`. Scalar extraction is limited to
the validated numerical formulas that evaluate the inertial-lift and Dean-drag
laws.

## Ownership decision

Aequitas owns the SI `Force` dimension and coherent `Newton` unit. CFDrs owns
the margination law and its equilibrium result. A consumer-local force wrapper
would duplicate the provider contract, while keeping the result as `f64` would
discard units at the public boundary.

## Verification

The Aequitas dimension-law test proves pressure-area composition into Newtons.
The cfd-1d and cfd-validation cell-separation tests preserve the typed
equilibrium position, force balance, direct formula inputs, Fahraeus/CFL/
rheology outputs, plasma-skimming partitioning, cross-junction geometry and
flow behavior, straight-versus-curved Dean drag, three-population behavior,
and model-geometry value semantics.
