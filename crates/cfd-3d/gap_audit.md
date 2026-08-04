# CFD-3D Gap Audit

## Aequitas and Eunomia audit (2026-07-29, CFDRS-AEQ-3D-01)

The public Venturi/FEM contract now constructs shared fluid properties through
`ConstantPropertyFluid` and Aequitas quantities. The cfd-3d audit found no
additional missing physical metric dimensions in the public solver result;
velocity, pressure, density, viscosity, heat capacity, conductivity, and
derived flow diagnostics are typed at their existing domain boundaries.

The FEM path is real-valued because its boundary validation and constitutive
ordering require Eunomia `RealField`. Complex Eunomia values remain applicable
only to phasor/spectral formula boundaries in other CFDrs consumers; no
imaginary-unit Aequitas quantity is needed for this real-valued flow contract.

The solver residuals found during this audit are closed: constrained DOFs are
removed before Krylov iteration and restored exactly, the configured relative
tolerance is used on every Picard linear system, P1 PSPG stabilization is
assembled once per element, and the fallback chain uses provider-backed
component factors with symbolic reuse. The Venturi blood validation remains
the value-semantic regression gate.

The historical API compatibility findings were revalidated against the
current peer tree with `cargo check -p cfd-3d --all-targets`; that gate is
green, so no cfd-3d compilation gaps remain in this audit.
