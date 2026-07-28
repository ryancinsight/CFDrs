# Cavitation-VOF physical metrics

## Scope

`cfd-3d` VOF cavitation configuration, bubble-dynamics construction, the
public cavitation time-step entry point, and their in-tree callers use
Aequitas quantities for values with a declared SI dimension.

## Decision

The public contracts carry:

- `SurfaceTension<f64>` for VOF surface tension and bubble surface tension;
- `Length<f64>` for bubble radius and mesh spacing at bubble-solver setup;
- `NumberDensity<f64>` for nuclei population density;
- `Time<f64>` for relaxation time, bubble updates, and cavitation steps;
- `Pressure<f64>` for vapor, ambient, and local bubble pressures;
- `MassDensity<f64>` for liquid and vapor densities;
- `Velocity<f64>` for sound speed.

Formula kernels still use the existing scalar Rayleigh-Plesset, damage, VOF,
and dense-field representations. Conversion is explicit at those boundaries,
so the dimensional contract remains visible at construction and call sites
without adding a second scalar facade.

Inception markers, void fractions, cavitation fractions, damage values, and
other normalized/model-specific statistics remain scalar. They do not have a
single SI dimension that Aequitas could represent without changing their
model contract.

## Rejected alternative

Retaining raw public fields with unit comments or adding typed accessors beside
the raw fields would preserve dimensional ambiguity and create two owners for
the same metric. The migration therefore changes the public Rust contracts and
updates every in-tree caller in one slice.

## Verification

The command-line local-provider overlay compiles `cfd-3d` test targets. The
focused cfd-3d Nextest run covers `cavitation_solver_validation`,
`robustness_tests`, and `vof_tests`: 83/83 passed. Source residue scans find
no raw typed VOF configuration assignments in the affected cfd-3d scope.

The ordinary standalone child build remains affected by duplicate local/git
provider identities and peer-owned lockfile dirt; those are integration
conditions recorded in `gap_audit.md`, not evidence against this metric
contract.
