# Microvascular blood Aequitas metrics

`CFDRS-AEQ-MET-33` closes the fixed-dimension storage boundary in
`cfd-core::physics::fluid::blood::FahraeuasLindqvist`.

## Typed public contract

The calculator stores:

- `Length<T>` for vessel diameter;
- `Dimensionless<T>` for feed and tube hematocrit;
- `DynamicViscosity<T>` for plasma and apparent viscosity.

The Pries and Secomb correlations operate on base scalars only at the
empirical formula boundary. Diameter is converted to micrometres there because
the published correlations use that unit. The cfd-1d adapter and PyO3 wrapper
construct the typed core values and extract scalars only at their existing
formula or external serialization boundaries.

## Eunomia compatibility

The calculator requires Eunomia `RealField` because the correlations compare
diameter and evaluate real powers and exponentials. Eunomia `Complex<T>` is not
an alternative representation for these ordered material properties, and no
imaginary-unit SI quantity is introduced. Complex values remain appropriate at
phasor, Bessel, Womersley, and spectral formula/storage boundaries.

## Verification

The exact increment passes cfd-core, cfd-1d, and cfd-python test-target checks,
focused blood-model Nextest (3/3), warning-denied cfd-core Clippy, cfd-core
doctests, no-default-features cfd-core rustdoc, targeted rustfmt, diff checks,
and the typed public-field residue scan.
