use crate::domain::model::NetworkBlueprint;
use cfd_core::error::Result;

/// Port for materializing a network blueprint into an output graph form.
pub trait GraphSink {
    /// Output graph representation.
    type Output;

    /// Build the output graph from a validated blueprint.
    fn build(&self, blueprint: &NetworkBlueprint) -> Result<Self::Output>;
}
