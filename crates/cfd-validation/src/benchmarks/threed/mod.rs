//! 3D benchmark submodules for CFD validation

pub mod bifurcation;
pub mod venturi;
pub mod serpentine;

pub use bifurcation::BifurcationFlow3D;
pub use venturi::VenturiFlow3D;
pub use serpentine::SerpentineFlow3D;
