//! Millifluidic channel construction.
//!
//! Generates watertight mesh geometry for microfluidic/millifluidic channels,
//! substrates, and junctions. Adapted from blue2mesh's extrusion pipeline
//! but using indexed mesh storage.

pub mod junction;
pub mod path;
pub mod profile;
pub mod substrate;
pub mod sweep;

pub use junction::JunctionType;
pub use path::ChannelPath;
pub use profile::ChannelProfile;
pub use substrate::SubstrateBuilder;
pub use sweep::SweepMesher;

pub mod branching;
pub mod serpentine;
pub mod venturi;

pub use branching::BranchingMeshBuilder;
pub use serpentine::SerpentineMeshBuilder;
pub use venturi::VenturiMeshBuilder;
