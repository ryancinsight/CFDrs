# ADR 0005: Type resistance model flow conditions

- Status: Proposed
- Class: [major]
- Date: 2026-11-21

Use Aequitas velocity, volumetric flow, reciprocal time, thermodynamic
temperature, pressure, and dimensionless Reynolds number in the flow
record. Pressure is absolute at the material query. A prescribed flow
leaves velocity absent so the model can divide by its effective area.
The original scalar formulae and material trait extract canonical SI
values at their existing numerical boundaries.

The Mach check uses the dimensionless ratio of speed magnitude and
sound speed. A test checks both signs at the 0.3 boundary, the adjacent
representable velocities above it, and rejection of a nonfinite velocity.
