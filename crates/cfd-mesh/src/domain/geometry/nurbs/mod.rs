//! # NURBS Curves and Surfaces
//!
//! Non-Uniform Rational B-Spline geometry, first-class in `cfd-mesh`.
//!
//! ## Modules
//!
//! | Module | Contents |
//! |---|---|
//! | [`knot`] | Knot vector validation, span search, span multiplicity |
//! | [`basis`] | Cox–de Boor recursive B-spline basis evaluation |
//! | [`curve`] | `NurbsCurve` — 3-D NURBS curve in homogeneous coordinates |
//! | [`surface`] | `NurbsSurface` — tensor-product NURBS surface |
//! | [`tessellate`] | Curvature-adaptive tessellation into triangle soup |
//!
//! ## Mathematical Background
//!
//! ### B-Spline Basis Functions — Cox–de Boor Recursion
//!
//! Given a knot vector `Ξ = {ξ₀, …, ξₘ}` and degree `p`, the i-th basis
//! function of degree `p` is:
//!
//! ```text
//! N_{i,0}(ξ) = 1  if  ξᵢ ≤ ξ < ξᵢ₊₁
//!            = 0  otherwise
//!
//! N_{i,p}(ξ) =  (ξ − ξᵢ)          N_{i,p−1}(ξ)
//!              ────────────────── · ──────────────
//!              (ξᵢ₊ₚ − ξᵢ)
//!
//!            +  (ξᵢ₊ₚ₊₁ − ξ)       N_{i+1,p−1}(ξ)
//!              ────────────────── · ───────────────
//!              (ξᵢ₊ₚ₊₁ − ξᵢ₊₁)
//! ```
//!
//! (Zero divided by zero is defined as 0.)
//!
//! ### Theorem — Partition of Unity
//!
//! For any knot span `[ξᵢ, ξᵢ₊₁)`:
//!
//! ```text
//! Σᵢ Nᵢ,ₚ(ξ) = 1    for all ξ
//! ```
//!
//! This guarantees that a B-spline curve (or surface) is an *affine
//! combination* of its control points, so translations and rotations of the
//! control polygon map exactly to translations and rotations of the curve.
//!
//! ### NURBS Rational Extension
//!
//! A NURBS curve with weights `wᵢ` is:
//!
//! ```text
//!        Σᵢ Nᵢ,ₚ(ξ) · wᵢ · Pᵢ
//! C(ξ) = ──────────────────────
//!          Σᵢ Nᵢ,ₚ(ξ) · wᵢ
//! ```
//!
//! ## Curvature-Adaptive Tessellation
//!
//! A span is subdivided only if the **angle between adjacent span normals**
//! exceeds a threshold `θ_max`.  This avoids unnecessary triangles in flat
//! regions and concentrates them at highly-curved areas — critical for smooth
//! CFD boundary interfaces.

pub mod knot;
pub mod basis;
pub mod curve;
pub mod surface;
pub mod tessellate;

pub use curve::NurbsCurve;
pub use surface::NurbsSurface;
pub use tessellate::TessellationOptions;
