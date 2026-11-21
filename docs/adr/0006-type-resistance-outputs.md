# ADR 0006: Type hydraulic resistance model outputs

- Status: Proposed
- Class: [major]
- Date: 2026-11-21

## Decision

Return Aequitas hydraulic resistance `R` and quadratic hydraulic
resistance `k` from the model and calculator. Preserve typed coefficients
through blueprint and network-edge assignments.

The signed pressure drop is `R Q + k Q |Q|`. With volumetric flow
`[m³/s]`, its terms have units `[Pa·s/m³][m³/s]` and
`[Pa·s²/m⁶][m³/s]²`. Multiplying `k` by `|Q|` gives effective
linear resistance where that is the existing model contract.

## Alternatives and verification

Returning a scalar pair leaves swapped coefficients and dimensionally
invalid addition undetectable. A typed pair keeps existing destructuring
and generic scalar precision. Geometry and the underlying correlation
evaluation preserve their existing numeric conventions.

Compile-fail tests reject mixing coefficients. The turbulent Darcy model
regression checks signed pressure at negative flow and verifies agreement
between the coefficient law and the model's effective resistance.
