//! Millifluidic channel construction.
//!
//! Generates watertight mesh geometry for microfluidic/millifluidic channels,
//! substrates, and junctions. Adapted from blue2mesh's extrusion pipeline
//! but using indexed mesh storage.

pub mod profile;
pub mod path;
pub mod sweep;
pub mod junction;
pub mod substrate;

pub use profile::ChannelProfile;
pub use path::ChannelPath;
pub use sweep::SweepMesher;
pub use junction::JunctionType;
pub use substrate::SubstrateBuilder;
