# cfd-core shared fluid properties

## Contract

`FluidProperties<T>` and `ConstantPropertyFluid<T>` are the shared storage
contracts for real-valued CFD material properties. Their public physical
fields use Aequitas quantities:

| Property | Contract type |
| --- | --- |
| Density | `MassDensity<T>` |
| Dynamic viscosity | `DynamicViscosity<T>` |
| Specific heat capacity | `SpecificHeatCapacity<T>` |
| Thermal conductivity | `ThermalConductivity<T>` |
| Speed of sound | `Velocity<T>` |

`PropertyBounds<T>` uses the matching quantity for every lower and upper
bound. `FluidProperties` derives `KinematicViscosity`, `ThermalDiffusivity`,
and `Dimensionless` Reynolds, Prandtl, Peclet, and Mach metrics. The
thermophysical subset delegates validation and diffusivity to Proteus through
its typed quantity constructor.

## Scalar boundaries

Base scalars are extracted only where the existing numerical contract requires
them: analytical formulas, dense field or mesh storage, FEM/GPU parameter
blocks, and explicit serialization. Callers construct the typed quantities at
the public boundary; no raw-field compatibility facade remains.

## Eunomia compatibility

The shared physical-property contracts require Eunomia `RealField` because
validation depends on ordering, positivity, and finite non-negative
thermophysical values. Eunomia complex scalars are therefore not widened into
these real material quantities. Complex values remain valid at phasor,
Womersley, Bessel, and spectral formula/storage boundaries; adding an
imaginary-unit SI property here would conflate complex response data with a
real material state.

## Remaining model boundary

`IdealGas`, temperature-dependent, non-Newtonian, and blood model structs still
store constitutive coefficients as raw generic scalars. Those values are model
parameters rather than the shared constant-property state and are the next
provider-first migration boundary. Their formula kernels must be converted
with the same explicit scalar-boundary rule before this residual is closed.

## Verification

The closure is tracked by `CFDRS-AEQ-MET-30` in `backlog.md` and `gap_audit.md`.
The focused value tests cover each derived metric and invalid-property error
contract; package checks and dependent test-target checks cover all migrated
constructors and callers.
