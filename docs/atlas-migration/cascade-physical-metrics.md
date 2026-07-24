# Cascade physical metric boundary

## Context

`cfd-3d::cascade` already used Aequitas internally to derive inlet velocity,
but discarded the types at the public configuration and result boundary. The
same API therefore admitted metres, pascals, cubic metres per second, and
metres per second as indistinguishable `f64` values.

## Decision

`CascadeChannelSpec` stores `Length` and `VolumetricFlowRate`; `CascadeConfig3D`
stores `Pressure`; and cascade results store `Pressure` and `Velocity`. The
existing JSON field names remain stable through explicit serde representation
types that materialize canonical SI scalars. Mesh coordinates, FEM boundary
conditions, and extracted scalar diagnostics convert once at their numerical
kernel boundary.

The generic venturi, bifurcation, and FEM internals remain scalar because they
operate on validated numerical buffers and are not part of this public cascade
DTO boundary. They are not duplicated or wrapped by a local compatibility
layer.

## Verification contract

- Typed constructors preserve canonical SI values at every migrated field.
- Existing cascade, adversarial, and cross-fidelity tests retain their
  pressure/shear/velocity value oracles after explicit extraction.
- Serde representation types preserve the established scalar wire shape.
- The focused package gate is pending because the current peer provider graph
  fails before `cfd-3d` source compilation: the local Coeus dependency resolves
  to `D:\atlas\repos\coeus\coeus-core\Cargo.toml`, while the checked-out
  manifest is at `D:\atlas\repos\coeus\crates\coeus-core\Cargo.toml`.
