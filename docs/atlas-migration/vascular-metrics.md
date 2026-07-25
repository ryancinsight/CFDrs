# Aequitas vascular metrics

## Context

The vascular modules exposed physical geometry, material properties, boundary
conditions, and derived results as unconstrained scalar values. That allowed a
radius, pressure, viscosity, flow, or time value to cross module boundaries
without a unit-bearing contract. The gap covered Womersley pulsatile flow,
vessel bifurcation networks, Murray geometry, and Olufsen structured-tree
resistance.

## Decision

Use Aequitas quantities at public physical boundaries and extract base scalars
only where an existing analytical kernel evaluates a dimensionless formula or
an algorithm-specific numerical recurrence.

- Womersley stores `Length`, `ReciprocalTime`, `MassDensity`, and
  `DynamicViscosity`; profile and flow results use `Velocity`, `Pressure`,
  `VolumetricFlowRate`, `PressureGradient`, and `HydraulicResistance`.
- Bifurcation vessels and networks store typed geometry, material properties,
  pressure, flow, resistance, inertance, compliance, area, and velocity.
- Murray bifurcation geometry stores `Length`, `Angle`, and
  `VolumetricFlowRate`; deviation, area ratio, and conservation error are
  `Dimensionless`.
- Olufsen stores the terminal radius as `Length` and returns structured-tree
  impedance as `HydraulicResistance`.
- PyO3 and validation wrappers convert at their boundary and retain their
  existing scalar-facing external contracts.

The Aequitas provider owns the named vascular dimensions. CFDrs owns only the
vascular equations and converts to base scalars at those equation boundaries.

## Rejected alternative

Retaining raw scalar fields with unit comments or adding parallel typed
accessors would preserve dimensional ambiguity and duplicate ownership. A
consumer-owned compatibility wrapper would also keep the old contract alive.

## Verification

- Aequitas provider: locked check and Nextest 28/28 passed in `446eb9f`.
- CFDrs `cfd-1d`: locked check and Nextest 731/731 passed with 3 skipped.
- CFDrs adversarial Womersley tests: 2/2 passed.
- CFDrs `cfd-python` and `cfd-validation`: locked checks pass.
- Remaining warning-denied Clippy failure is the pre-existing
  `cfd-math::matrix_zeros` dead-code warning; `leto-ops` reports existing
  dead helpers under ordinary warning output.
