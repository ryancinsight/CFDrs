//! 3D Geometry submodules for CFD validation

pub mod bifurcation;
pub mod serpentine;
pub mod venturi;

pub use bifurcation::Bifurcation3D;
pub use serpentine::Serpentine3D;
pub use venturi::Venturi3D;
