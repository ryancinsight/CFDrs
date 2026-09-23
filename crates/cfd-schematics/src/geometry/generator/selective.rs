mod build;
mod builder;
mod geometry;
mod path_geometry;
mod pending;
mod primitive;
mod request;
mod routing;

pub use primitive::{
    PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
    create_primitive_selective_tree_geometry, create_primitive_selective_tree_geometry_from_spec,
};

pub use build::create_selective_tree_geometry;
pub use request::{CenterSerpentinePathSpec, SelectiveTreeRequest, SelectiveTreeTopology};

// Private aliases: `builder/*` and `primitive/annotation` reach these two DTOs
// as `super::…`, which is where they sat before the split.
use geometry::SelectiveTreeGeometry;
use pending::PendingVenturiPath;

#[cfg(test)]
mod tests;
