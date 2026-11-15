# ADR 0001: Advance the Eunomia and Moirai provider contracts

Status: Accepted  \
Date: 2026-11-15  \
Change class: [arch] [major]  \
Board item: [ATLAS-VERSION-COHERENCE-SWEEP-2026-09-30](../../backlog.md)

## Context

CFDrs binds its public scalar contract to Eunomia's `RealField`
(`cfd-core::scalar`, `cfd-io`'s CSV types) and resolves its parallel
dispatch through Moirai's runtime. The stack's provider cohort has moved
past both lines:

- Eunomia 0.12.0 is the current provider head. Rust treats the
  `RealField` trait of Eunomia 0.11 and the same-named trait of Eunomia
  0.12 as distinct nominal contracts, so a generic consumer whose bound
  still names the 0.11 trait no longer satisfies CFDrs's public bounds
  even when its scalar type implements both versions.
- Moirai 0.12.0 removes the required `StateWord::get_mut` method. CFDrs
  implements no `StateWord` and calls no removed method, so no CFDrs
  source changes with the move.

The requirement also gates the rest of the stack. The providers' current
default heads require Moirai 0.12 (`apollo-fft` 0.29.5 and
`coeus-core` 0.12.0 name `moirai-runtime ^0.12.0`), so every consumer
that resolves those providers with a `^0.11.0` requirement cannot
re-resolve onto the current provider heads at all.

Two alternatives were considered:

1. Stay on Eunomia 0.11 and Moirai 0.11 and keep the current lock. The
   provider cohort then splits across two Eunomia identities and two
   Moirai identities inside one workspace, and the consumers that
   resolve CFDrs inherit the same split.
2. Require `>=0.11.0, <0.13.0` for Eunomia and let the lock choose. The
   range admits both provider lines, but the lock selects the newest
   compatible revision anyway, so the re-exported and bound traits move
   exactly as they do under the selected option. The wider range only
   hides the move from the requirement.

## Decision

Move the Eunomia requirement to 0.12.0, the Moirai runtime requirement to
0.12.0, and widen the Leto and Leto-ops requirements to
`>=0.46.4, <0.49.0` so the 0.48 line resolves. Regenerate the lock
outside the Atlas dependency overlay, so every provider row names one
revision per provider repository and the `source = "git+..."` rows
survive.

No CFDrs source changes with the move. The removed provider method is
unused in this repository.

## Consequences and migration

- Code generic over the Eunomia 0.11 `RealField` bound no longer
  satisfies CFDrs's public type bounds. A consumer updates its Eunomia
  bound and imports to 0.12 and selects the same provider cohort before
  it uses CFDrs's public generic types. Consumers already using concrete
  scalar types keep the same CFDrs API after their dependency graph
  selects the new cohort.
- The public contract change is a [major] change class; the release
  tooling assigns the version at release.
- No compatibility trait, forwarding wrapper, or duplicate scalar
  contract is added: two Eunomia identities under one name would keep
  the split the provider owns.
