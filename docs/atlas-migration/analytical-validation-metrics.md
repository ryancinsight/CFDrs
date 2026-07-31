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

The same current-head scan found raw fixed-dimension fields in the Couette,
Poiseuille, Stokes, Taylor-Green, Blasius, and non-Newtonian analytical
configurations, plus the legacy `analytical_benchmarks` module. Those remain
separate dependency-ordered metric items. Dense sampled fields, coordinates,
interpolation tables, dimensionless coefficients, and runtime-exponent
formula parameters remain explicit scalar boundaries.

## Verification

The Womersley public-field residue scan, targeted Rustfmt, and `git diff
--check` pass. The pinned cfd-validation test-target check is currently
blocked before this crate by peer-dirty `cfd-math` unresolved `leto-ops`
symbols in `linear_solver/block_preconditioner.rs:26-28` and the missing
`factor_sparse_with_symbolic` method at line 769.
