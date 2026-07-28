# Cavitation physical metrics

## Scope

`cfd-core` cavitation primitives, the standalone `cfd-3d` closure seam, and
`cfd-3d` VOF cavitation configuration, bubble-dynamics construction, the public
cavitation time-step entry point, and their in-tree callers use Aequitas
quantities for values with a declared SI dimension.

## Decision

The public contracts carry:

- `SurfaceTension<f64>` for VOF surface tension and bubble surface tension;
- `Length<f64>` for bubble radius and mesh spacing at bubble-solver setup;
- `NumberDensity<f64>` for nuclei population density;
- `Time<f64>` for relaxation time, bubble updates, and cavitation steps;
- `Pressure<f64>` for vapor, ambient, and local bubble pressures;
- `MassDensity<f64>` for liquid and vapor densities;
- `Velocity<f64>` for sound speed.
- `DynamicViscosity<f64>`, `ThermodynamicTemperature<f64>`, and `Energy<f64>`
  for Rayleigh-Plesset parameters and sonoluminescence estimates;
- `Dimensionless<f64>` for cavitation numbers and normalized model fractions;
- `Frequency<f64>` for regime-classifier acoustic frequency;
- `Angle<f64>` and `Volume<f64>` for Venturi geometry and cavity outputs;
- `ThermalDiffusivity<f64>` for nuclei transport diffusion;
- `MassDensityRate<f64>` for phase-transfer rates.
- typed water surface-tension, vapor-pressure, and nuclei-density constants
  in the cfd-core cavitation constants module.

Formula kernels still use scalar arithmetic at equation, dense-field, GPU, or
serialization boundaries. Conversion is explicit at those boundaries, so the
dimensional contract remains visible at construction and call sites without
adding a second scalar facade.

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

The command-line local-provider overlay compiles cfd-core and cfd-3d test
targets. cfd-core Nextest passes 202/202 with no skips. The broad cfd-3d
Nextest run passes 291/292; every cavitation, VOF, closure, robustness, and
validation test passes, while the sole timeout is the pre-existing
`test_venturi_blood_flow` runtime residual at 30.663 seconds against the
committed 30-second budget. The runtime residual is tracked under
`CFDRS-RUNTIME-001`, not this metric contract. The prior MET-24 focused
evidence remains 83/83. Source residue scans find no raw typed VOF
configuration assignments or placeholder collapse-rate implementation in the
affected cfd-3d scope.

The ordinary standalone child build remains affected by duplicate local/git
provider identities and peer-owned lockfile dirt; those are integration
conditions recorded in `gap_audit.md`, not evidence against this metric
contract.
