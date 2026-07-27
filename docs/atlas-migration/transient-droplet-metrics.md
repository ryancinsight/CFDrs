# Transient droplet physical metrics

## Contract

The public droplet-tracking boundary carries physical and normalized values as
Aequitas quantities:

- `DropletInjection::volume`, `DropletSnapshot::total_volume`, and split-policy
  minimum child volume use `Volume`.
- `DropletInjection::injection_time` and `DropletTrackingState::time` use
  `Time`.
- Injection, boundary, position, and occupancy coordinates use `Dimensionless`
  because they are normalized channel fractions.
- Flow-weighted policy thresholds use `Dimensionless`.

Quantity values remain in canonical SI base units. The simulator converts to
the scalar representation only in advection, finite-length occupancy, and
flow-weighted splitting formulas. Internal branch state carries the same typed
volume and normalized-position values so the public snapshot cannot diverge
from its conserved state.

The composition engine's independent event and time-state contracts are not
part of this slice; their remaining raw fields stay tracked as the next
transient-composition audit item rather than being silently relabeled closed.

## Rejected alternative

Keeping droplet volume and time as `T` and relying on SI comments was rejected:
the compiler could not distinguish a volume from a time or normalized position
at injection and snapshot call sites. A compatibility constructor was rejected
because the public contract must migrate atomically.

## Verification

- `cargo check -p cfd-1d --tests --offline` passed before the peer cfd-math
  graph shifted.
- Direct rustfmt and targeted `git diff --check` pass.
- The focused Nextest command is currently blocked by the peer-owned
  `cfd-math::BlockDiagonalPreconditioner` and `SimplePreconditioner`
  `apply_to` return-type mismatch (`LetoError` versus `cfd_core::error::Error`).
