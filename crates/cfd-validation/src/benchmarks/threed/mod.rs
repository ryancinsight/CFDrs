//! 3D benchmark submodules for CFD validation

pub mod bifurcation;
pub mod serpentine;
pub mod venturi;

pub use bifurcation::BifurcationFlow3D;
pub use serpentine::SerpentineFlow3D;
pub use venturi::VenturiFlow3D;
