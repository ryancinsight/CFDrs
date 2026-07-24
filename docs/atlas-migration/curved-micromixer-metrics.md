# Curved-channel and micromixer geometry metrics

## Context

`cfd-1d` documented the `ChannelType::Curved` radius and micromixer hydraulic
diameter/path length in metres but stored them as raw generic or `f64`
scalars. The values therefore crossed public construction and storage
boundaries without a dimension-level contract.

## Decision

Use Aequitas `Length` for:

- `ChannelType::Curved::radius`;
- `Micromixer::hydraulic_diameter`; and
- `Micromixer::length`.

Resistance formulas extract base scalars at their numerical boundary. The
string-keyed `Micromixer::set_parameter` API remains a scalar configuration
boundary and immediately reconstructs a typed `Length`; it does not create a
second stored geometry representation.

## Rejected alternative

Retaining scalar fields with metre comments would preserve dimensional
ambiguity. Adding typed accessors beside the scalar fields would create two
owners for the same geometry. Both alternatives are rejected.

## Verification contract

The locked cfd-1d test check and configured Nextest cover curved geometry,
micromixer resistance, constructor validation, dynamic parameter updates,
component validation, and adversarial callers. The acceptance result is
731/731 tests passed with 3 configured skips. Doctest, Rustdoc, Clippy, and
runtime-budget evidence remains in `gap_audit.md`.
