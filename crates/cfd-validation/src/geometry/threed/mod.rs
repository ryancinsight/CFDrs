//! 3D Geometry submodules for CFD validation

pub mod bifurcation;
pub mod venturi;
pub mod serpentine;

pub use bifurcation::Bifurcation3D;
pub use venturi::Venturi3D;
pub use serpentine::Serpentine3D;
