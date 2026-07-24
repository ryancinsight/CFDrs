# Component geometry and volume metrics

## Context

The public `cfd-1d` component API documented metres and square metres but
stored channel, membrane, and organ geometry as raw scalars. The shared
`Component::volume` contract also returned an untyped scalar, and
`ChannelProperties` repeated the same untyped channel dimensions.

## Decision

Use Aequitas quantities at the component boundary:

- channel, membrane, and organ lengths, widths, heights, diameters,
  roughnesses, and pore radii use `Length`;
- component area methods return `Area`;
- the shared component volume method returns `Volume`;
- `ChannelProperties` stores `Length` values.

Resistance models still receive base scalars at their numerical formula
boundary. The string-keyed factory and `set_parameter` methods are explicit
scalar configuration boundaries and immediately construct or update typed
quantities. No scalar compatibility constructors or fields remain.

## Rejected alternative

Keeping the component fields scalar would leave dimensional assignment errors
possible in the primary public model. Adding parallel typed accessors while
retaining raw fields would create two owners for the same geometry and is not
retained.

## Verification contract

The component unit, property, adversarial, and membrane-parity tests preserve
resistance, monotonicity, and volume value semantics. The locked `cfd-1d`
check, Nextest, doctests, Rustdoc, and targeted residue scan are the acceptance
gates. Existing provider-graph warnings, lint failures, and documentation link
warnings remain recorded in `gap_audit.md` rather than being hidden.
