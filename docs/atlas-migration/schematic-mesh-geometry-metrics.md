# Schematic mesh geometry metrics

## Context

`cfd-schematic-mesh` exposed runtime mesh configuration and emitted layout
geometry as unit-suffixed scalars. Angles, chip and cavity dimensions, wall
clearances, CSG overlap fractions, centerline coordinates, and diameters could
therefore cross the public mesh boundary without a dimensional type check.

## Decision

Use Aequitas quantities for the runtime mesh contract:

- `PipelineConfig` carries `Angle`, `Length`, and `Dimensionless` values;
- `ShellPipelineConfig` carries cavity height, chip height, and cavity mid-plane
  as `Length` values;
- `SegmentCenterline` carries coordinates and diameter as `Length` values.

The mesh provider still consumes millimetre/radian scalars at the explicit
trigonometry, routing, primitive, and CSG formula boundaries. The serialized
`cfd-schematics` interchange payload remains a separate unit-labelled format
boundary and is not treated as a second runtime owner of these values.

No complex or imaginary physical quantity applies: mesh geometry is real under
Eunomia's real scalar contract. Eunomia complex values remain valid only for
genuine phasor/Fourier data elsewhere in the stack.

## Rejected alternative

Keeping raw public fields would preserve dimensional ambiguity. Adding typed
accessors beside the scalar fields would create two owners and retain the
untyped path. Both alternatives are rejected.

## Verification contract

The focused `cfd-schematic-mesh` package check and value-semantic Nextest are
required. The standalone provider graph is currently blocked before package
compilation because the locked Moirai revision requests a `mnemosyne` package
that is not available in the local git checkout; this external provider
residual is recorded in `gap_audit.md` and does not change the geometry
contract.
