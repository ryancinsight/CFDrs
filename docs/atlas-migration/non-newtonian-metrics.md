# Non-Newtonian Aequitas metrics

`CFDRS-AEQ-MET-32` closes the fixed-dimension storage boundary in
`cfd-core::physics::fluid::non_newtonian`.

## Typed public contract

The Bingham, Casson, Carreau–Yasuda, Power-law, and Herschel–Bulkley models
carry fixed-dimension parameters with Aequitas:

- `MassDensity` for density;
- `Pressure` for yield stress;
- `DynamicViscosity` for Newtonian viscosity parameters;
- `ReciprocalTime` for reference shear rate;
- `Time` for Carreau relaxation time;
- `Dimensionless` for flow and Yasuda exponents;
- `ThermodynamicTemperature` for reference temperature;
- `MolarEnergy` for activation energy;
- `SpecificHeatCapacity`, `ThermalConductivity`, and `Velocity` for thermal
  and acoustic material properties.

Consistency indices remain scalar formula data. In `K · γ̇^n`, their unit is
`Pa·s^n`, which changes with the runtime exponent; assigning one fixed
Aequitas dimension would be incorrect. Apparent-viscosity methods therefore
extract the typed shear-rate and coefficient values only at the formula
boundary. The legacy `Fluid` and `NonNewtonianFluid` trait scalar interfaces
remain separate compatibility boundaries for numerical solver inputs.

## Eunomia compatibility

These constitutive models require Eunomia `RealField` for ordered validation and
real-valued powers and exponentials. Eunomia `Complex<T>` and an imaginary-unit
SI material quantity are not applicable to density, stress, viscosity, or
thermal-property state. Complex values remain at phasor, Bessel, Womersley, and
spectral formula/storage boundaries.

## Verification

The exact increment passes cfd-core test-target check, focused non-Newtonian
Nextest (4/4), warning-denied all-targets Clippy, cfd-core doctests (3/3),
no-default-features cfd-core rustdoc, targeted rustfmt, diff checks, and the
typed public-field residue scan. Default-feature rustdoc remains blocked by an
unrelated peer `hephaestus-wgpu` `Send + Sync` bound error at
`application/elementwise_seam.rs:413`.
