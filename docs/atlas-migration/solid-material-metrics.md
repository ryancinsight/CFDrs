# cfd-core solid material metrics

## Contract

`cfd-core::physics::material::SolidProperties` and `ElasticSolid` are public
material contracts. Their unit-bearing values are represented by Aequitas:

| Metric | Contract type |
| --- | --- |
| Density | `MassDensity<T>` |
| Young's and shear modulus | `Pressure<T>` |
| Poisson's ratio | `Dimensionless<T>` |
| Thermal conductivity | `ThermalConductivity<T>` |
| Specific heat capacity | `SpecificHeatCapacity<T>` |
| Thermal expansion | `ReciprocalTemperature<T>` |

The isotropic relation

`G = E / (2 * (1 + nu))`

is evaluated through typed Aequitas operations and returns `Pressure<T>`.
There is no raw scalar mirror of the public fields or trait methods.

## Eunomia boundary

The trait remains constrained to `T: FloatElement + Copy`, matching the
real-valued isotropic material model. Eunomia `Complex<T>` is therefore not a
valid scalar for this contract, and no imaginary-unit quantity is introduced.
Complex constitutive material data would require a separate provider-backed
contract with an explicit phasor semantic; it must not be represented by a
real metric or an untyped complex scalar.

## Migration

Callers constructing `ElasticSolid` wrap canonical SI values with
`Quantity::from_base`. Numerical consumers use `into_base()` only where a
formula explicitly requires a scalar; this material contract's shear-modulus
formula remains typed end-to-end.

## Verification

The acceptance evidence is tracked by `CFDRS-AEQ-MET-28` in `backlog.md` and
the corresponding entry in `gap_audit.md`.
