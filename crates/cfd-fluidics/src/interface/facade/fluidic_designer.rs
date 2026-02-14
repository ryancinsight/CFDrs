use crate::application::use_cases::NetworkGenerationService;
use crate::domain::model::NetworkBlueprint;
use crate::infrastructure::adapters::Cfd1dGraphSink;
use cfd_1d::NetworkGraph;
use cfd_core::error::Result;

pub struct FluidicDesigner {
    generator: NetworkGenerationService<Cfd1dGraphSink>,
}

impl FluidicDesigner {
    #[must_use]
    pub fn new() -> Self {
        Self {
            generator: NetworkGenerationService::new(Cfd1dGraphSink),
        }
    }

    pub fn generate(&self, blueprint: &NetworkBlueprint) -> Result<NetworkGraph<f64>> {
        self.generator.generate(blueprint)
    }
}

impl Default for FluidicDesigner {
    fn default() -> Self {
        Self::new()
    }
}
