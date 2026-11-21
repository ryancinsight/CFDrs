# ADR 0003: Type blood reference constants

- Status: Proposed
- Class: [major]
- Date: 2026-11-21

The blood reference module stores empirical viscosity, yield stress, density,
time, thermal properties, sound speed, hematocrit, and vessel diameter as
Aequitas quantities. Body temperature carries an absolute temperature and
the Andrade coefficient carries a temperature difference. Their product
with the reciprocal-temperature difference is dimensionless.

Model constructors cast tabulated values once into their declared scalar
precision. Raw solvers and Python boundaries extract canonical SI values
explicitly. The migration keeps the empirical values, model equations, and
numeric encodings. The obsolete rounded Casson square-root coefficient has
no caller. The solver derives its viscosity term from the typed
infinite-shear viscosity.

The density compile-fail example rejects a viscosity assignment. The
Fåhræus–Lindqvist test checks the cutoff and adjacent floating-point
diameters.

## Shared whole-blood density

Cavitation loading, Venturi optimization, and cross-fidelity validation
read the same `MassDensity` constant from the blood reference module.
The optimizer exports that constant directly and computes dynamic pressure
as `ρv²/2` using dimensional arithmetic. Lower-level scalar protocols read
the canonical SI value at their boundaries. The book demonstrations and
network benchmark use the same density. Cell-specific reference
densities and independent numerical fixtures retain their own values.
