# Aequitas analytical-validation metrics

## Context

The public `cfd-validation::analytical::WomersleyFlow` configuration stored
radius, density, viscosity, angular frequency, and pressure-gradient amplitude
as raw generic scalars even though the canonical `cfd-1d` Womersley evaluator
already accepted typed values. The boundary therefore discarded dimensions
before validation and direct reuse.

## Decision

Store the fixed-dimension Womersley configuration with Aequitas:

- `Length` for radius;
- `MassDensity` for density;
- `DynamicViscosity` for viscosity;
- `ReciprocalTime` for angular frequency; and
- `PressureGradient` for the driving gradient.

Return `Dimensionless`, `Velocity`, `Length`, `Pressure`, and
`VolumetricFlowRate` from the corresponding derived metrics. Extract base
scalars only for the analytical Bessel/formula implementation and the
`AnalyticalSolution` mesh-coordinate boundary. The canonical `cfd-1d`
Womersley profile remains the single evaluator; no consumer wrapper or
parallel scalar API is retained.

## Eunomia compatibility

The public model remains constrained by Eunomia `RealField` because its
validation and classification metrics require ordered real values. The
canonical evaluator may use complex Bessel/phasor intermediates internally,
but those values do not cross the Aequitas physical boundary. Eunomia's
`Complex<T>` is the real-plus-quadrature representation; it does not define a
separate imaginary-unit physical quantity, and this model does not need one.

## Follow-on audit

The same current-head scan found raw fixed-dimension fields in the Blasius and
non-Newtonian analytical configurations, plus the legacy
`analytical_benchmarks` module. Those remain separate dependency-ordered metric
items. Dense sampled fields, coordinates, interpolation tables, dimensionless
coefficients, and runtime-exponent formula parameters remain explicit scalar
boundaries.

## Verification

The Womersley public-field residue scan, targeted Rustfmt, and `git diff
--check` pass. The pinned cfd-validation test-target check is currently
blocked before this crate by peer-dirty `cfd-math` unresolved `leto-ops`
symbols in `linear_solver/block_preconditioner.rs:26-28` and the missing
`factor_sparse_with_symbolic` method at line 769.

## Couette and Poiseuille

`CouetteFlow` stores wall velocity, gap height, pressure gradient, and
dynamic viscosity as Aequitas `Velocity`, `Length`, `PressureGradient`, and
`DynamicViscosity`. Its shear rate, wall stress, and Reynolds number return
`ReciprocalTime`, `Pressure`, and `Dimensionless` respectively.

`PoiseuilleFlow` stores its maximum velocity, characteristic length, pressure
gradient, and viscosity with the same typed boundary. Its velocity and
Reynolds calculations return `Velocity` and `Dimensionless`. The runtime
geometry is closed by `PoiseuilleFlowRate`: plates return
`AreaPerTime` (`m²/s`, flow per unit width), while pipes return
`VolumetricFlowRate` (`m³/s`). This preserves the distinction instead of
coercing both formulas to one scalar dimension.

The formula and `AnalyticalSolution` mesh-coordinate boundaries extract base
values explicitly. These validation contracts require Eunomia `RealField`;
complex intermediates and quadrature values are not physical-unit inputs for
this family, so no imaginary-unit quantity is introduced.

## Stokes flow around a sphere

`StokesFlow` stores sphere radius, free-stream velocity, dynamic viscosity,
and density as Aequitas `Length`, `Velocity`, `DynamicViscosity`, and
`MassDensity`. Its drag law returns `Force`, the drag coefficient and Reynolds
number return `Dimensionless`, and the spherical stream function returns
`VolumetricFlowRate` because its definition is velocity times area.

The analytical vector and pressure methods retain raw mesh-coordinate values
at the `AnalyticalSolution` boundary. The physical configuration and direct
metrics remain real-valued under Eunomia `RealField`; no complex or
imaginary-unit quantity is introduced.

## Taylor-Green vortex

`TaylorGreenVortex` stores its characteristic length, velocity, kinematic
viscosity, and density with Aequitas `Length`, `Velocity`,
`KinematicViscosity`, and `MassDensity`. Its Reynolds number, decay rate, and
enstrophy return `Dimensionless`, `ReciprocalTime`, and
`ReciprocalTimeSquared`.

The explicit `TaylorGreenDimension` branch prevents a dimensional mismatch in
kinetic energy. Two-dimensional solutions return
`TaylorGreenKineticEnergy::PerDepth(Force)` because their energy is per unit
depth; three-dimensional solutions return
`TaylorGreenKineticEnergy::Volumetric(Energy)`. The analytical and mesh
boundaries extract base scalars only for trigonometric formulas and coordinate
arrays. The 3D DNS benchmark converts its typed analytic result at the
benchmark/report boundary.

The model remains real-valued under Eunomia `RealField`; no complex or
imaginary-unit quantity applies. Aequitas owns `ReciprocalTimeSquared` rather
than a consumer-local unit alias.
