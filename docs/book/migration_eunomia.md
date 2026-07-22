# Eunomia: Numeric Trait Unification

Every Atlas crate depends on [`eunomia::RealField`] for scalar arithmetic and
on [`eunomia::ComplexField`] when complex numbers appear.  Eunomia is the
**trait frontier** — a one-crate abstraction that the rest of the stack builds
on, replacing `num-traits`, `num-complex`, and the various `num-derive`
utilities that historically scatter across the dependency graph.

## The RealField Trait

`eunomia::RealField` is the minimum surface a type must implement to participate
in Atlas:

```rust
pub trait RealField:
    Copy
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + PartialOrd
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(x: f64) -> Self;
    fn to_f64(self) -> f64;
    fn sqrt(self) -> Self;
    fn abs(self) -> Self;
    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
}
```

The trait is **deliberately small**.  Anything beyond the minimum (`sin`,
`exp`, `powf`, `ln`, ...) lives in extension traits, so a type only pays for
the surface it actually uses.

## Companion Traits

```rust
pub trait FloatElement: RealField + ... { /* f32/f64 only rules */ }
pub trait IntElement: ... { /* i32/i64 only rules */ }
pub trait ComplexField: ... { /* z = re + i*im rules */ }
```

These traits **carve the trait frontier in two**: `FloatElement` for floating-point
work, `IntElement` for indexing arithmetic, `ComplexField` for spectral
methods.

## ZST and PhantomData Bounds

Atlas crates use ZSTs and `PhantomData` to encode capability requirements
**at the type level**, with zero runtime cost:

```rust
pub struct Solver<F: FloatElement, B: ComputeBackend> {
    inner: PhantomData<(F, B)>,
    state: Vec<F>,
}
```

`PhantomData<(F, B)>` is **zero-sized** — it compiles away to nothing — but
forces the compiler to verify that `Solver<f32, Cpu>` and `Solver<f64, Wgpu>`
are distinct types.  This is how Atlas delivers **monomorphization without
runtime dispatch**: the same source code generates specialized machine code
per `(F, B)` pair, with no virtual call sites.

## Migration From num-traits

| Legacy | Atlas |
|---|---|
| `T: num_traits::Float` | `T: eunomia::FloatElement` |
| `num_complex::Complex64` | `eunomia::ComplexField` impl |
| `from_f64` via `as` cast | `RealField::from_f64` |
| `fn add(a: f64, b: f64)` | RealField bound at the API point |

CFDrs normally ports its solver entry points like this:

```rust
// legacy: fn step(&mut self, dt: f64) — flat f64 only
// atlas:
pub fn step<F: FloatElement>(&mut self, dt: F) -> Result<(), CfdError> {
    self.state.iter_mut().for_each(|x| *x = *x + dt);
    Ok(())
}
```

The function now compiles once per `F` with no dynamic dispatch and no `as`
cast at the boundary.

## Further Reading

- [`eunomia` source](../../../eunomia/crates/)
- [Atlas Dependency Map](appendix_dependencies.md)
