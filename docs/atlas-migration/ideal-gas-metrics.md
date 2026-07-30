# Ideal-gas Aequitas metrics

`CFDRS-AEQ-MET-34` closes the fixed-dimension metric boundary in
`cfd-core::physics::fluid::newtonian::IdealGas`.

## Typed public contract

`IdealGas` stores its constitutive parameters with Aequitas:

- `SpecificHeatCapacity` for the specific gas constant, constant-pressure
  heat capacity, and conductivity coefficient;
- `DynamicViscosity` for the reference viscosity;
- `ThermodynamicTemperature` for the reference temperature; and
- `TemperatureDifference` for the Sutherland offset.

The `Fluid` trait retains its scalar condition signature as a separate legacy
boundary. `IdealGas::properties_at` converts those inputs immediately to
`ThermodynamicTemperature` and `Pressure`, then returns the existing typed
`FluidState` metrics.

The gas constant and conductivity coefficient intentionally use
`SpecificHeatCapacity`: both have the coherent SI dimension `J/(kg·K)`, and
Aequitas does not need a semantic alias for values with the same physical
dimension.

## Formula boundaries

Base scalar extraction is limited to the ideal-gas density equation,
Sutherland viscosity, kinetic-theory conductivity, and ideal-gas sound-speed
formula. The stored model and returned fluid state remain typed.

## Eunomia compatibility

The model requires Eunomia `RealField` because state validity uses ordering and
positive temperatures/pressures, and Sutherland's law and sound speed are
real-valued constitutive formulas. Eunomia `Complex<T>` and an imaginary-unit
Aequitas material quantity are not applicable. Complex values remain at
phasor, Bessel, Womersley, or spectral boundaries.

## Verification

The exact increment passes the cfd-core test-target check, selected newtonian
Nextest 16/16, warning-denied Clippy, cfd-core doctests 3/3,
no-default-features rustdoc, targeted rustfmt, staged diff checks, and the
typed public-field residue scan. Default-feature rustdoc remains subject to
the unrelated peer `hephaestus-wgpu` `Send + Sync` error at
`application/elementwise_seam.rs:413`.
