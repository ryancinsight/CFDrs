# Hemolysis exposure metric boundary

## Context

`cfd-1d` exposed wall shear stress and exposure duration as untyped `f64`
arguments to its Giersiepen and Taskin models. `HemolysisExposure` repeated
the same untyped pair in public fields. The resulting hemolysis index and
cavitation potential are dimensionless model values and do not need a
physical-quantity wrapper.

## Decision

Use Aequitas `Pressure` for shear stress and `Time` for exposure duration in
the `cfd-1d` public functions and `HemolysisExposure` carrier. Keep the
Giersiepen formula constants and Taskin formula in their existing `cfd-core`
or local model owners. Convert typed values to canonical SI scalars only at
the `cfd-core::HemolysisModel` call and at scalar analysis/reporting
boundaries.

The migration updates all in-repository callers, examples, and validation
tests. No scalar forwarding overload is retained.

## Verification contract

- Zero, negative, reference-value, monotonicity, and model-divergence tests
  retain their analytical hemolysis-index oracles.
- Flow analysis and optimization reporting convert existing scalar storage
  into typed boundary values before invoking the models.
- Package format, checks, Nextest, warning-denied Clippy, doctests, and Rustdoc
  are rerun after the peer Coeus manifest path is repaired; the current Cargo
  graph stops before CFDrs source compilation at
  `D:\atlas\repos\coeus\coeus-core\Cargo.toml`.
