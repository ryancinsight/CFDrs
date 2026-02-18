use crate::application::ports::GraphSink;
use crate::domain::model::NetworkBlueprint;
use crate::domain::rules::BlueprintValidator;
use cfd_core::error::Result;

pub struct NetworkGenerationService<S: GraphSink> {
    sink: S,
}

impl<S: GraphSink> NetworkGenerationService<S> {
    #[must_use]
    pub fn new(sink: S) -> Self {
        Self { sink }
    }

    pub fn generate(&self, blueprint: &NetworkBlueprint) -> Result<S::Output> {
        BlueprintValidator::validate(blueprint)?;
        self.sink.build(blueprint)
    }
}
