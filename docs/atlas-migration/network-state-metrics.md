# Network state physical metrics

## Decision

The public `cfd-1d` network state carries the Aequitas dimensions already
owned by the provider:

- `Network::pressures` and pressure accessors use `Pressure<T>`.
- `Network::flow_rates` and flow accessors use `VolumetricFlowRate<T>`.
- `NetworkState::time` uses `Time<T>`.
- Network analysis results use the corresponding pressure, flow, velocity,
  reciprocal-time, power, resistance, residence-time, and dimensionless
  quantities.

Scalar extraction is restricted to matrix assembly, resistance and transport
formulas, mesh/GPU-facing calculations, and explicit serialized/reporting
boundaries. In-tree callers construct quantities at setter boundaries; no raw
scalar compatibility facade remains.

## Residual classification

`Network::residuals` and `PrimarySolveDiagnostics::last_residual_norm` remain
scalar solver diagnostics. Their unit is determined by the assembled equation
and its scaling, and the current solver admits both dimensional and
normalized residual systems. Assigning one Aequitas physical dimension would
misrepresent at least one valid equation configuration. The values are not
public physical measurements and are therefore outside this metric slice.

## Verification contract

The acceptance oracle is value-semantic: typed state construction and
conversion preserve the previous base values, solver flow convergence uses
base values only inside its norm calculation, and analysis formulas produce
the same scalar results after typed wrapping. The focused locked check,
Nextest, doctest, Rustdoc, and warning-denied Clippy results are recorded in
`gap_audit.md` when the increment closes.
