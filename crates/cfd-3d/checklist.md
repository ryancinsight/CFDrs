# CFD-3D Update Checklist

## CFDRS-AEQ-MET-44 typed turbulence metrics — closed 2026-07-31

- [x] Migrate the shared turbulence contract to Aequitas
      `KinematicViscosity<T>` and `SpecificEnergy<T>` results.
- [x] Migrate every canonical turbulence closure and its value-semantic
      regression assertions; keep scalar extraction at formula boundaries.
- [x] Record the Eunomia boundary: real turbulence uses real scalars;
      complex values remain for genuine phasor/Fourier fields and no imaginary
      physical unit is introduced.
- [x] Record the exact local locked-gate blocker: the shared Atlas overlay
      requests a provider lock refresh before rustc starts.

## Historical compilation-error reconciliation — closed 2026-07-31

- [x] Audit the historical compilation-error list against the current source.
- [x] Replace the obsolete `num_traits::Float`/`ComplexField` ambiguity item
      with the current Eunomia `RealField`/`FloatElement` contracts.
- [x] Confirm the `IndexedMesh` vertex/face calls, serpentine builder, dense
      matrix bounds, and iterator scalar conversions compile in the current
      cfd-3d source.
- [x] Verify the public turbulence path resolves to the canonical
      `physics::turbulence` implementations rather than no-op placeholders.
- [x] Record the current verification boundary: the pinned audit records
      `cargo check -p cfd-3d --all-targets` green and the current source passes
      metadata/diff checks; this checkout cannot start a locked compile while
      the shared Atlas overlay requests a lockfile refresh for its mutable
      provider patches.
