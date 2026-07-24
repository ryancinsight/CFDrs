# Channel geometry and network edge metrics

## Context

The public `cfd-1d` network model documented metres, square metres, and
hydraulic diameters but represented channel cross-sections, channel lengths,
edge areas, and edge-property geometry as raw scalars. The same values then
crossed resistance, transient transport, analysis, and example boundaries
without a type-level dimension check.

## Decision

Use Aequitas quantities for the public channel/network geometry contract:

- `CrossSection` dimensions, custom area, and hydraulic diameter use `Length`
  and `Area`;
- `ChannelGeometry.length` uses `Length`;
- `Edge.area` uses `Area`;
- `EdgeProperties.length`, `area`, and `hydraulic_diameter` use `Length` and
  `Area`.

Resistance, junction-loss, transient-transport, and analysis formulas extract
base scalars only at their existing numerical kernels. Blueprint conversion,
dynamic property construction, examples, and test fixtures are explicit
scalar-to-quantity boundaries. No parallel scalar fields or compatibility
constructors remain.

## Rejected alternative

Keeping the geometry fields scalar would preserve dimensional ambiguity. Adding
typed accessors beside the raw fields would create two owners for the same
geometry and leave downstream callers free to use the untyped path. Both
alternatives are rejected.

## Verification contract

The locked `cfd-1d` check and configured Nextest cover the geometry kernels,
network assembly, transient transport, analyzers, blueprint conversion,
examples, property tests, and adversarial tests. The acceptance result is
729/729 tests passed with 3 configured skips. Doctest, Rustdoc, Clippy, and
runtime-budget results remain recorded in `gap_audit.md` with their exact
limits; no provider-graph or lint warning is treated as a false green.
