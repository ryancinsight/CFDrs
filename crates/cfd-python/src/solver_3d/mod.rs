//! 3D CFD solver `PyO3` wrappers

mod bifurcation;
mod poiseuille;
mod venturi;

pub use bifurcation::*;
pub use poiseuille::*;
pub use venturi::*;
