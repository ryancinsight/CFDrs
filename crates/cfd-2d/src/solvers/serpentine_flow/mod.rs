//! 2D serpentine channel flow solver with mixing efficiency validation
//!
//! This module implements validated solvers for serpentine (sinuous) channels,
//! which are critical for microfluidic mixing and sample preparation.
//!
//! # Physics Background
//!
//! ## Serpentine Channel Design
//!
//! Serpentine channels enhance mixing through:
//! 1. **Advection**: Fluid streams are stretched and folded
//! 2. **Secondary flow**: Dean vortices in curves increase transverse mixing
//! 3. **Residence time**: Longer path length increases reaction time
//!
//! ## Mixing Mechanisms
//!
//! In laminar microfluidic flows (Re << 1), mixing is diffusion-limited:
//!
//! **Diffusion time scale:**
//! ```text
//! t_diff = w² / (4D)  [seconds, from Fick's law]
//! ```
//!
//! where:
//! - w = channel width \[m]
//! - D = diffusion coefficient [m²/s]
//!
//! **Advection time scale:**
//! ```text
//! t_adv = L / u  [seconds, from advection]
//! ```
//!
//! **Mixing length:**
//! ```text
//! L_mix = u · t_diff = u · w² / (4D)
//! ```
//!
//! For complete mixing across channel width: L_adv > L_mix
//!
//! ## Peclet Number
//!
//! Ratio of advection to diffusion:
//! ```text
//! Pe = u·w / D
//! ```
//!
//! - Pe << 1: Diffusion-dominated mixing
//! - Pe >> 1: Mixing limited by diffusion length
//!
//! ## Mixing Efficiency
//!
//! Quantified by intensity of segregation (ISO):
//! ```text
//! ISO = ∫(c - c_mean)² dV / ∫c_max² dV
//! ```
//!
//! where:
//! - ISO = 1 at inlet (completely unmixed)
//! - ISO = 0 at outlet (perfectly mixed)
//!
//! # Validation Strategy
//!
//! 1. **Advection-diffusion**: Solve transport equation with known inlet conditions
//! 2. **Richardson extrapolation**: Grid convergence study
//! 3. **Literature**: Compare against experimental mixing data
//! 4. **Analytical**: Use advection-diffusion analytical solutions
//!
//! # References
//!
//! - Hardt, S. & Schönfeld, F. (2003). "Microfluidic technologies for miniaturized
//!   analysis systems". SpringerLink
//! - Squires, T.M. & Quake, S.R. (2005). "Microfluidics: Fluid physics at the
//!   nanoliter scale". Reviews of Modern Physics, 77(3), 977
//! - Yao, Z., et al. (2014). "Numerical study of mixing in microchannels with
//!   oscillatory flow". Microfluidics and Nanofluidics, 16(1), 145-155
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

mod geometry;
mod mixing;
mod solver;
mod validation;

pub use geometry::SerpentineGeometry;
pub use mixing::{AdvectionDiffusionMixing, SerpentineMixingSolution};
pub use solver::SerpentineSolver2D;
pub use validation::{SerpentineValidationResult, SerpentineValidator};
