# Network metadata physical metrics

## Contract

`cfd-1d` network metadata carries declared physical quantities through the
public domain boundary. Node pressure uses Aequitas `Pressure`, node
temperature uses `ThermodynamicTemperature`, network volume uses `Volume`, and
pressure and temperature ranges use pairs of their corresponding quantities.
The quantity values are stored in Aequitas canonical SI base units.

The free-form `NodeProperties::metadata` map remains scalar because its keys
do not declare a physical dimension. Treating that map as one named quantity
would make unknown metadata appear dimensionally meaningful and would be a
more damaging contract than retaining the explicitly untyped extension point.

## Boundary rule

Typed quantities remain at public constructors, fields, and builder setters.
Scalar extraction is permitted only at numerical formula, mesh/GPU, or
explicit serialization boundaries. This slice has no compatibility wrapper and
does not retain parallel raw scalar fields.

## Rejected alternative

Keeping `T` fields and documenting their SI interpretation was rejected because
the compiler could not distinguish pressure, temperature, and volume at call
sites. Adding a local dimension enum for the free-form metadata map was also
rejected because it would invent semantics not present in its keys.

## Verification

- `cargo check -p cfd-1d --offline`
- `cargo nextest run -p cfd-1d --offline`: 735 passed, 3 skipped
- `cargo test --doc -p cfd-1d --offline`: 8 passed, 3 ignored
- `cargo clippy -p cfd-1d --all-targets --offline -- -D warnings`: reaches the
  peer-owned `cfd-math` bridge and stops on its missing-docs warning; no
  diagnostic originates in this slice.
