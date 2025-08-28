//! Mesh domain module with clean separation of concerns

mod connectivity;
mod operations;
mod quality;
mod statistics;
mod types;

pub use connectivity::{Connectivity, EdgeConnectivity, FaceConnectivity};
pub use operations::{MeshOperations, MeshTransform};
pub use quality::{MeshQuality, QualityMetrics};
pub use statistics::MeshStatistics;
pub use types::{Element, Mesh, MeshMetadata};

use crate::Result;
use nalgebra::RealField;
