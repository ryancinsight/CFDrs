//! 2D CFD solver `PyO3` wrappers.
//!
//! This module exposes 2D finite volume solvers and 1D resistance models
//! for validation against Python CFD packages.
//!
//! ## Sub-modules
//!
//! | Module | Solvers |
//! |--------|---------|
//! | [`poiseuille`] | `Poiseuille2DSolver`, `Poiseuille2DResult` |
//! | [`venturi`] | `VenturiSolver2D/1D`, `VenturiResult2D/1D` |
//! | [`cavity`] | `CavitySolver2D`, `CavityResult2D` |
//! | [`bifurcation`] | `BifurcationSolver2D`, `BifurcationResult2D`, `TrifurcationSolver2D`, `TrifurcationResult2D` |
//! | [`serpentine`] | `SerpentineSolver1D`, `SerpentineResult1D` |

mod bifurcation;
mod cavity;
mod poiseuille;
mod serpentine;
mod venturi;

pub use bifurcation::*;
pub use cavity::*;
pub use poiseuille::*;
pub use serpentine::*;
pub use venturi::*;
