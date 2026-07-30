# Temperature-model Aequitas metrics

`CFDRS-AEQ-MET-31` closes the temperature-dependent model storage boundary in
`cfd-core::physics::fluid::temperature`.

## Typed public contract

`PolynomialViscosity`, `ArrheniusViscosity`, `AndradeViscosity`, and
`SutherlandViscosity` store representable physical parameters with Aequitas:

- `MassDensity` for density;
- `ThermodynamicTemperature` for absolute/reference evaluation temperatures;
- `TemperatureDifference` for Arrhenius, Andrade, and Sutherland offsets;
- `ReciprocalTemperature` for the polynomial density expansion coefficient;
- `DynamicViscosity`, `SpecificHeatCapacity`, `ThermalConductivity`, and
  `Velocity` for the corresponding constitutive properties.

The calculation methods accept `ThermodynamicTemperature`. The legacy scalar
condition in the shared `Fluid` trait is converted once at the implementation
boundary because that trait is a separate migration surface.

Polynomial viscosity coefficients remain `Vec<T>` formula data. The coefficient
at order `i` has units of `Pa·s/K^i`, so one quantity type cannot represent the
whole vector without falsely claiming a common dimension. Scalar extraction is
therefore explicit and limited to polynomial evaluation and the Proteus
temperature-response formula.

## Eunomia compatibility

These models require Eunomia `RealField` for ordering, positivity checks, and
real-valued constitutive functions. Eunomia `Complex<T>` and an imaginary-unit
Aequitas quantity are not applicable: complex values belong at phasor or
spectral formula/storage boundaries, not to density, viscosity, or thermal
property state.

## Verification

The exact increment passes cfd-core test-target check, warning-denied
Clippy/all-targets, temperature-model Nextest (7/7), cfd-core doctests (3/3),
targeted rustfmt, diff checks, and the typed public-field residue scan. The
remaining model-property work is limited to ideal-gas, non-Newtonian, and
blood constitutive coefficients.
