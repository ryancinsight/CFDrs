# ADR 0002: CFDrs release scope for 1.0

Status: Accepted  \
Date: 2026-11-15  \
Change class: [arch]  \
Board item: CFDRS-RELEASE-LADDER

## Context

CFDrs is the stack's fluid-momentum and mass owner (1-D networks, 2-D/3-D
finite-volume, GPU Poisson/advection/diffusion via Hephaestus, validation
against manufactured and benchmark oracles). The board carries seven open
delivery items (JFNK recovery, three SIMPLEC performance items, the Python
surface pair, one Apollo-blocked example gate), an in-progress lint-floor
increment, and a dozen gap-audit hygiene items; the structural-improvement
run is fully landed. No milestone ladder orders them toward a stable 1.0,
and no record states what 1.0 freezes.

## Decision

Adopt a seven-rung ladder to 1.0. The ladder opens `backlog.md`; every open
item carries exactly one `milestone:` assignment.

- **M-newton**: bounded Newton-Krylov recovery lands
  (ATLAS-CFDRS-JFNK-RECOVERY-110, P0). Exit: the JFNK acceptance oracle passes
  at head; no open P0 solver item.
- **M-solver-perf**: problem-scaled SIMPLEC targets hold
  (ATLAS-CFDRS-RUNTIME-109, ATLAS-CFDRS-SOLID-PRESSURE-CACHE-108,
  ATLAS-CFDRS-PRESSURE-CACHE-102). Exit: each item's benchmark oracle
  passes without repeated matrix clones or unscaled termination.
- **M-validation**: the validation tier is honest (CFDRS-CAVITY-VALIDATION
  review, CFDRS-GA-006, CFDRS-GA-003). Exit: the cavity suite asserts
  falling centerline error under refinement against Ghia Table I,
  convergence studies assert observed order, and no Python domain logic
  remains outside `crates/cfd-python/tests`.
- **M-python**: the Python surface is typed and publishable
  (CFDRS-GA-014A, ATLAS-PUBLISH-001-CFDRS-PYPI; CFDRS-GA-014 delivered as
  exit evidence). Exit: GIL released around Rust-only work, complete typed
  stub with PEP 561 marker, value-semantic installed-wheel tests, and
  trusted-publishing wheels through the shared Atlas caller.
- **M-blocked-external**: the Apollo-blocked example gate
  (ATLAS-CFDRS-APOLLO-PLAN-SCRATCH-104) stays a recorded hold outside the
  blocking path; it reopens when the provider export lands.
- **M-hygiene**: gates and records are mechanized (CFDRS-LINT-FLOOR-001,
  CFDRS-GA-002, GA-005, GA-007, GA-008, GA-009, GA-010, GA-011, GA-015,
  GA-018). Exit: lint counts strictly decrease, root examples compile
  under a gate, `expect("expected value")` count is zero, figures verify
  content, one board/checklist/gap file each, book chapters numbered and
  non-stub, `missing_docs` denied per crate, and the PR scope filter runs
  fewer jobs on source-only changes.
- **M-1.0**: stable fluids release (**authority-gated**). Exit: `cargo
  semver-checks`, publish dry-run, `cargo doc`, the full workspace
  nextest/doctest suite green, and zero open correctness items.
  **Veto point:** release authority sets the API-freeze contents; this
  milestone never forces `publish = true` without it.
- **M-ecosystem**: consumer-owned coupling (prometheus source terms via
  harmonia, kwavers acoustic coupling). Exit: consumer receipts on their
  own boards; does not gate 1.0.

### Veto points

1. The 1.0 API-freeze contents (release authority sets them).
2. The nalgebra transitional boundary: the board's residuals record
   nalgebra-sparse storage behind Leto facades; full removal is its own
   item, not a hidden 1.0 condition.

## Rejected alternatives

- Freezing 1.0 while a P0 solver item is open: JFNK recovery is the
  correctness backstop for stiff networks; a stability promise cannot
  cover its absence.
- Counting the Apollo-blocked gate as a 1.0 blocker: the provider export
  is outside this repo's control; it stays a recorded hold.

## Overturning evidence

A consumer failure traceable to a frozen API, or a validation oracle
failing at head, reopens the affected rung.
