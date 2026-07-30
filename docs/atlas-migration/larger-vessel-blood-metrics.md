# Larger-vessel blood metrics

## Scope

CFDRS-AEQ-MET-35 moves the fixed-dimension state of the cfd-core Casson,
Carreau–Yasuda, and Cross blood models onto Aequitas quantities. The change
covers the model constructors, public fields, constitutive accessors, direct
tests, cfd-2d and cfd-3d examples, cfd-validation formulas, and cfd-python
serialization accessors.

The typed fields are:

- MassDensity
- Pressure
- DynamicViscosity
- Dimensionless
- Time
- ReciprocalTime
- SpecificHeatCapacity
- ThermalConductivity
- Velocity

Fluid and NonNewtonianFluid results remain typed as DynamicViscosity. The
apparent-viscosity formula continues to accept the formula's scalar shear-rate
representation and returns its scalar formula value; the trait and property
boundaries carry the result as an Aequitas quantity. Temperature correction
now accepts a typed absolute ThermodynamicTemperature.

## Boundary policy

Base-scalar extraction is limited to constitutive formulas, mesh or solver
adapters, and explicit Python serialization. Fixed-dimension model state does
not retain parallel scalar accessors or compatibility wrappers. Consistency
indices whose dimensions depend on a runtime flow exponent remain formula
parameters rather than being assigned a false fixed unit.

The blood-rheology models use Eunomia RealField because validation, ordering,
positivity, and real-valued empirical powers are part of their contract.
Eunomia Complex<T> and an imaginary-unit SI material property do not apply to
this ordered rheology state. Complex values remain at phasor, spectral,
Bessel, and Womersley formula or storage boundaries where the governing
quantity is genuinely complex.

## Verification

The focused cfd-core check, blood-model Nextest, warning-denied Clippy,
doctests, targeted rustfmt, typed-field residue scan, and staged diff checks
are the acceptance gates for this slice. Consumer package checks are reported
with their exact provider or peer diagnostic when the shared Atlas workspace
prevents compilation before the changed consumer is reached.
