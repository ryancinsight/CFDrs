//! 1D bifurcation and trifurcation solvers for microfluidic networks
//!
//! This module implements validated solvers for branching flow networks with:
//! - Bifurcation/trifurcation junction conservation equations
//! - Non-Newtonian blood fluid models (Casson, Carreau-Yasuda)
//! - Validated against analytical solutions (Huo & Kassab, Fung)
//! - Complete documentation in code for reproducibility
//!
//! # Physics Background
//!
//! ## Bifurcation Junction Equations
//!
//! At a bifurcation junction, conservation laws apply:
//!
//! **Mass Conservation (continuity):**
//! ```text
//! Q_parent = Q_1 + Q_2
//! ```
//!
//! **Pressure Distribution (dynamic junction):**
//! For incompressible flow, the pressure at the junction satisfies:
//! ```text
//! P_parent = P_1 + f_1(Q) = P_2 + f_2(Q)
//! ```
//! where f_i is the pressure drop in each daughter branch.
//!
//! ## Blood Flow Specifics
//!
//! Blood exhibits non-Newtonian behavior crucial in bifurcations:
//! - Shear rate varies across the junction
//! - Viscosity is shear-rate dependent (Casson, Carreau-Yasuda models)
//! - Fåhræus-Lindqvist effect in small vessels (D < 300 μm)
//!
//! # Validation Strategy
//!
//! Simulations are validated using:
//! 1. **Conservation checks**: Mass conservation error must be < 1e-10
//! 2. **Analytical comparisons**: Against Poiseuille-based solutions
//! 3. **Convergence studies**: Spatial convergence with Richardson extrapolation
//! 4. **Literature benchmarks**: Huo & Kassab (2012) branching models
//!
//! # References
//!
//! - Huo, Y., & Kassab, G. S. (2012). "Intraspecific scaling laws of vascular trees"
//! - Fung, Y. C. (1993). "Biomechanics: Mechanical Properties of Living Tissues"
//! - Zamir, M. (1992). "The Physics of Pulsatile Flow"

pub mod junction;
pub mod network_solver;
pub mod validation;

pub use junction::{BifurcationJunction, BifurcationSolution, TrifurcationJunction, TrifurcationSolution};
pub use network_solver::{BifurcationConfig, BifurcationNetworkSolver};
pub use validation::{BifurcationValidationResult, BifurcationValidator};
