# ADR 0003 · Type blood reference constants

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
