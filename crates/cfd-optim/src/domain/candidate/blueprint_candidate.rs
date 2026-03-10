use crate::domain::OptimizationGoal;
use crate::error::OptimError;
use cfd_schematics::{
    BlueprintTopologyFactory, BlueprintTopologySpec, ChannelSpec, NetworkBlueprint, NodeSpec,
    TopologyOptimizationStage,
};
use serde::{Deserialize, Serialize};

use super::OperatingPoint;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintCandidate {
    pub id: String,
    pub blueprint: NetworkBlueprint,
    pub operating_point: OperatingPoint,
}

impl BlueprintCandidate {
    #[must_use]
    pub fn new(
        id: impl Into<String>,
        blueprint: NetworkBlueprint,
        operating_point: OperatingPoint,
    ) -> Self {
        Self {
            id: id.into(),
            blueprint,
            operating_point,
        }
    }

    pub fn from_topology_spec(
        id: impl Into<String>,
        topology: &BlueprintTopologySpec,
        operating_point: OperatingPoint,
    ) -> Result<Self, OptimError> {
        let blueprint =
            BlueprintTopologyFactory::build(topology).map_err(OptimError::InvalidParameter)?;
        Ok(Self::new(id, blueprint, operating_point))
    }

    pub fn topology_spec(&self) -> Result<&cfd_schematics::BlueprintTopologySpec, OptimError> {
        self.blueprint.topology_spec().ok_or_else(|| {
            OptimError::InvalidParameter(format!(
                "blueprint candidate '{}' is missing topology metadata",
                self.id
            ))
        })
    }

    #[must_use]
    pub fn inferred_goal(&self) -> OptimizationGoal {
        match self
            .blueprint
            .lineage()
            .map(|lineage| lineage.current_stage)
            .unwrap_or(TopologyOptimizationStage::AsymmetricSplitResidenceSeparation)
        {
            TopologyOptimizationStage::AsymmetricSplitResidenceSeparation => {
                OptimizationGoal::AsymmetricSplitResidenceSeparation
            }
            TopologyOptimizationStage::AsymmetricSplitVenturiCavitationSelectivity => {
                OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity
            }
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement => {
                OptimizationGoal::InPlaceDeanSerpentineRefinement
            }
        }
    }

    #[must_use]
    pub fn treatment_channel_ids(&self) -> Vec<String> {
        self.blueprint.treatment_channel_ids()
    }

    #[must_use]
    pub fn blueprint(&self) -> &NetworkBlueprint {
        &self.blueprint
    }

    #[must_use]
    pub fn nodes(&self) -> &[NodeSpec] {
        &self.blueprint.nodes
    }

    #[must_use]
    pub fn channels(&self) -> &[ChannelSpec] {
        &self.blueprint.channels
    }

    #[must_use]
    pub fn node(&self, id: &str) -> Option<&NodeSpec> {
        self.blueprint
            .nodes
            .iter()
            .find(|node| node.id.as_str() == id)
    }

    #[must_use]
    pub fn channel(&self, id: &str) -> Option<&ChannelSpec> {
        self.blueprint
            .channels
            .iter()
            .find(|channel| channel.id.as_str() == id)
    }

    #[must_use]
    pub fn into_blueprint(self) -> NetworkBlueprint {
        self.blueprint
    }
}
