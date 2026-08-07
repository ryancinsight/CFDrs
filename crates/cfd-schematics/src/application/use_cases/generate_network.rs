use crate::application::ports::GraphSink;
use crate::domain::model::NetworkBlueprint;
use crate::domain::rules::BlueprintValidator;
use cfd_core::error::Result;

/// Use case that validates a blueprint and materializes it through a sink.
pub struct NetworkGenerationService<S: GraphSink> {
    sink: S,
}

impl<S: GraphSink> NetworkGenerationService<S> {
    /// Create a service bound to the given graph sink.
    #[must_use]
    pub fn new(sink: S) -> Self {
        Self { sink }
    }

    /// Validate the blueprint, then build its output graph.
    pub fn generate(&self, blueprint: &NetworkBlueprint) -> Result<S::Output> {
        BlueprintValidator::validate(blueprint)?;
        self.sink.build(blueprint)
    }
}
