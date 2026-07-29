# 2D cell-tracking physical metrics

## Decision

The public `cfd-2d::solvers::cell_tracking` contract carries physical values
as Aequitas quantities. The tracker accepts typed positions, velocities, time
steps, fluid properties, and channel geometry. Routing fractions and lift
coefficients remain `Dimensionless` quantities. Dense solver fields and the
analytical cell-force formulas continue to use base scalars only after an
explicit extraction at that numerical boundary.

The trajectory sample is a named `TrackedPosition` containing `Length` x/y
coordinates and `Time`, rather than a heterogeneous scalar array. This keeps
the time coordinate from being confused with a spatial coordinate and gives
callers a typed serialization contract.

## Boundary map

| Boundary | Typed contract | Scalar extraction |
| --- | --- | --- |
| Velocity-field interpolation | `Length` → `(Velocity, Velocity)` | Staggered-grid interpolation |
| Tracker configuration | `DynamicViscosity`, `MassDensity`, `Length`, `Velocity` | Drag, lift, buoyancy, and wall formulas |
| Cell identity/state | `Length`, `Velocity`, `MassDensity` | Particle integration state |
| Trajectory samples | `Length`, `Time` | Serialization only when a scalar wire format is selected |
| Bifurcation routing | `Length`, `Dimensionless` | Pries formula boundary |

The tracker is a real-valued particle model. Eunomia complex values do not
represent any public cell-tracking metric; no imaginary unit or complex
quantity extension is required. Complex-valued consumers must keep their
phasor or spectral values at their own formula/storage boundary.

## Verification oracle

The analytical Poiseuille and asymmetric-bifurcation tests remain the
behavioral oracle. They assert domain exit, trajectory sample production,
centerline stability, and the CTC-versus-RBC routing relation after the type
cutover. A residue scan covers the public cell-tracking fields and trait
surface; any remaining raw values are limited to the Pries formula DTO or
private numerical state.
