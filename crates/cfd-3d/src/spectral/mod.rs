//! Spectral methods for 3D CFD simulations.
//!
//! This module provides high-order spectral methods following
//! proper separation of concerns.
//!
//! ## Theorem — Exponential Convergence for Analytic Solutions
//!
//! **Theorem (Spectral Convergence, Gottlieb & Orszag 1977)**: For a function f
//! that is analytic in a strip |Im(z)| < δ around the real line, the N-mode
//! spectral approximation fₙ satisfies:
//!
//! ```text
//! ‖f - fₙ‖_∞ = O(e^{-δN})
//! ```
//!
//! where δ > 0 depends on the width of the analyticity strip. This contrasts
//! with algebraic O(N^{-p}) convergence of p-th order finite differences.
//!
//! ## Theorem — Parseval's Identity (Energy Conservation)
//!
//! **Theorem (Parseval 1799)**: For a periodic sequence {fⱼ}_{j=0}^{N-1} and
//! its DFT {F̂ₖ}:
//!
//! ```text
//! Σⱼ |fⱼ|² = (1/N) Σₖ |F̂ₖ|²
//! ```
//!
//! **Implication**: The transform is energy-preserving (unitary up to scaling).
//! Round-trip transform/inverse must recover energy to machine precision:
//! ‖f - IDFT(DFT(f))‖ < N·ε_mach.
//!
//! ## Theorem — Aliasing & Nyquist (Shannon 1949)
//!
//! **Theorem**: A signal with bandwidth B can be exactly reconstructed from
//! samples taken at rate f_s > 2B (Nyquist rate). For DFT of N samples:
//! - Modes k = 0, ..., N/2-1 are representable (Nyquist modes).
//! - A frequency k > N/2 aliases to frequency k - N.
//!
//! **Dealiasing Rule**: For convolution products of order p (e.g., quadratic
//! nonlinearity p=2), the minimum dealiased grid size is N_dealiased ≥ ⌈3N/2⌉
//! (the "3/2 rule"). This is enforced in the pseudospectral Navier-Stokes solver.
//!
//! ## Theorem — Gauss-Lobatto Orthogonality (Chebyshev)
//!
//! **Theorem**: The Chebyshev polynomials {Tₙ(x)}_{n=0}^N are orthogonal with
//! respect to the weight w(x) = 1/√(1-x²) on [-1,1]:
//!
//! ```text
//! ∫_{-1}^1 Tₙ(x) Tₘ(x) / √(1-x²) dx = cₘ π/2 · δ_{nm}
//! ```
//!
//! where c₀ = 2, cₙ = 1 for n ≥ 1. The Gauss-Lobatto quadrature at N+1 points
//! xⱼ = cos(jπ/N) integrates polynomials of degree ≤ 2N-1 exactly.
//!
//! ## Invariants
//! - Mode count N must satisfy N ≥ 2 (minimum for meaningful transform).
//! - For FFT efficiency, N should be a power-of-2 or have small prime factors.
//! - After inverse transform, max round-trip error must be < N·2.2e-16.
//! - Periodic pseudospectral DNS and seeded forcing use Apollo FFTs with a 2/3 de-aliasing filter.

pub mod basis;
pub mod chebyshev;
pub mod diagnostics;
pub mod dns;
pub mod forcing;
pub mod fourier;
pub mod poisson;
pub mod solver;

pub use basis::{BasisFunction, SpectralBasis};
pub use chebyshev::ChebyshevPolynomial;
pub use diagnostics::{
    enstrophy_spectrum, kinetic_energy_spectrum, probe_signal_spectrum, temporal_autocorrelation,
    EnstrophySpectrum, KineticEnergySpectrum, ProbeSignalSpectrum, TemporalAutocorrelation,
};
pub use dns::{PeriodicPseudospectralDns3D, PeriodicPseudospectralDnsConfig};
pub use forcing::{
    BandLimitedRandomPhaseForcing3D, BandLimitedRandomPhaseForcingConfig,
    TimeResampledBandLimitedForcing3D, TimeResampledBandLimitedForcingConfig,
};
pub use fourier::{FourierTransform, SpectralDerivative};
pub use poisson::{PoissonBoundaryCondition, PoissonSolver};
pub use solver::{SpectralConfig, SpectralSolution, SpectralSolver};
