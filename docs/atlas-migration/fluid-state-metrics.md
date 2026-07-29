# cfd-core FluidState metrics

## Contract

`cfd-core::physics::fluid::FluidState` carries the real-valued fluid state
through Aequitas quantities:

| Metric | Contract type |
| --- | --- |
| Density | `MassDensity<T>` |
| Dynamic viscosity | `DynamicViscosity<T>` |
| Specific heat capacity | `SpecificHeatCapacity<T>` |
| Thermal conductivity | `ThermalConductivity<T>` |
| Sound speed | `Velocity<T>` |
| Kinematic viscosity | `KinematicViscosity<T>` |
| Thermal diffusivity | `ThermalDiffusivity<T>` |
| Reynolds, Prandtl, Peclet, and Mach numbers | `Dimensionless<T>` |

The `Fluid`, `ConstantFluid`, `NonNewtonianFluid`, and `CompressibleFluid`
provider seams return the same typed quantities. Thermophysical diffusivity
delegates to Proteus through its typed `ThermophysicalProperties` constructor.

## Scalar boundaries

Scalar extraction is limited to formula inputs, dense solver or mesh storage,
FEM/GPU adapters, and explicit serialization. The cfd-1d and cfd-3d callers
convert at those boundaries with `into_base()`; they do not reconstruct raw
physical values beside the typed contracts.

The older `ConstantPropertyFluid` storage and `FluidProperties` constructor
family remain a separate raw-scalar compatibility boundary. They are not
represented as closed by this `FluidState` increment.

## Eunomia compatibility

This state model is real-valued and constrained by Eunomia `RealField` where
the equations require it. Eunomia `Complex<T>` is not a valid substitute for
these real physical state quantities. Complex values remain at their existing
Womersley, Bessel, and spectral formula or storage boundaries; introducing an
imaginary-unit SI quantity here would conflate phasor data with real fluid
properties.

## Solver behavior

The 2D channel adapter propagates field-solver errors and hemolysis validation
errors instead of substituting reference-trace metrics. Reverse directed flow
is solved as a magnitude under the field solver's geometric boundary contract,
then its sign is restored on the extracted flow metric. Pressure and resistance
remain field-derived values.

## Verification

The closure is tracked by `CFDRS-AEQ-MET-29` in `backlog.md` and
`gap_audit.md`. The exact acceptance gates are the cfd-core test-target check,
cfd-core Nextest, cfd-core doctests, warning-denied cfd-core and cfd-2d
Clippy, the cfd-2d Nextest suite, and dependent cfd-1d/cfd-3d/
cfd-validation test-target checks.
