# Proteus temperature-response ownership

## Context

`cfd-core::physics::fluid::PolynomialViscosity` historically evaluated the
linear thermal-expansion density law locally:

`rho(T) = rho_ref * (1 - beta * (T - T_ref))`.

Kwavers uses the same response family for tissue properties. Proteus now owns
the dimensional, validated, statically dispatched temperature-response law.
Retaining the CFDrs arithmetic would create a second constitutive-law owner.

## Decision

`cfd-core` delegates the density response to Proteus through its existing
thermophysical anti-corruption module. The CFD expansion coefficient `beta`
maps to the Proteus linear-response coefficient `-beta`; specific heat and
conductivity use `ConstantResponse` zero-sized strategies.

`PolynomialViscosity::calculate_density` becomes fallible because Proteus
rejects invalid reference temperatures, non-finite coefficients, invalid
thermophysical reference properties, and negative evaluated densities.
Callers propagate the existing `cfd_core::Error::InvalidInput` boundary.

The law remains monomorphized over `T`. Its GAT state borrows the evaluation
temperature, and its response set contains no allocation or dynamic dispatch.

## Rejected alternative

Keeping the scalar formula in CFDrs avoids a fallible public method but leaves
constitutive arithmetic duplicated and permits negative or non-finite density
states. A downstream wrapper around the scalar formula has the same ownership
defect.

## Verification

The focused tests compare the delegated result with an independent closed-form
linear oracle under a four-rounding native-precision bound and assert that a
coefficient producing negative density is rejected with Proteus property and
constraint diagnostics.
