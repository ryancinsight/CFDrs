# Workspace Hardening Audit — 2026-04-13

## Scope

This audit snapshot covers the executable policy and code paths touched in the first hardening pass:

- workspace nextest policy in `.config/nextest.toml`
- `crates/apollo/crates/apollofft`
- report/asset review contract at the workspace root

It is intentionally limited to verified findings gathered from the current tree. It does not claim complete coverage of every crate.

## Confirmed Findings

### 1. Test timeout policy conflicted with the required 30 second ceiling

Observed state:

- `.config/nextest.toml` contained overrides granting up to 120 seconds for selected tests.

Correction applied:

- removed all per-test timeout overrides
- retained a single default `slow-timeout = { period = "15s", terminate-after = 2 }`

Result:

- the workspace policy is now a strict 30 second nextest budget unless a future commit changes it explicitly

### 2. `apollofft` constructors admitted invalid dimensions through raw scalar APIs

Observed state:

- `FftPlan1D::new` and `FftPlan1D::with_precision` accepted raw `usize`
- `FftPlan2D` and `FftPlan3D` accepted raw dimensions directly
- validated domain types `Shape1D`, `Shape2D`, and `Shape3D` already existed but were not the public construction contract for plans

Correction applied:

- plan constructors now require validated shape types
- cache construction paths now accept validated shapes and build plans from them
- convenience helpers construct validated shapes before planning

Result:

- zero-sized and otherwise invalid plan shapes are rejected at the domain boundary instead of being admitted into the planning layer

### 3. `apollofft` 1D theorem/test coverage lagged behind the 3D implementation

Observed state:

- `plan1d.rs` had no local theorem statement or proof sketch
- 1D plan code lacked direct property coverage for linearity, Parseval, precision-profile ordering, and caller-owned buffer equivalence

Correction applied:

- added theorem/proof/assumptions/failure-mode rustdoc to `FftPlan1D`
- added deterministic and property-based tests for:
  - linearity
  - Parseval identity
  - roundtrip recovery
  - caller-owned buffer equivalence
  - mismatched buffer rejection
  - precision-profile error ordering

Result:

- 1D FFT behavior is now documented and executable as an invariant-bearing module rather than only an implementation detail

## Remaining Verified Gaps

These gaps were confirmed but not fully remediated in this pass:

- many non-FFT crates still contain long-running tests whose algorithms or fixture sizes need reduction rather than timeout exceptions
- theorem/proof coverage is still uneven outside `apollofft`
- report/image pipelines do not yet consistently emit review manifests
- visual inspection remains a manual release task and is not yet centrally recorded for all generated assets

## Immediate Next Targets

1. reduce or decompose the formerly long-running bifurcation, trifurcation, GA, and venturi tests until they pass the strict nextest policy without exceptions
2. propagate validated constructor policy into additional public APIs that still accept invalid raw dimensions or physical parameters
3. wire the report/figure manifest contract into `cfd-optim` reporting and example output generation
4. create a release-gating visual review ledger by opening generated SVG/PNG/PDF outputs and recording pass/fail observations
