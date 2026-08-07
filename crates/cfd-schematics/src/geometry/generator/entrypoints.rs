use super::super::metadata::{BlueprintRenderHints, MetadataContainer};
use super::super::types::SplitType;
use super::GeometryGenerator;
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::NetworkBlueprint;
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};
use aequitas::systems::si::quantities::Length;

/// Configuration for metadata generation.
#[derive(Debug, Clone, Default)]
pub struct MetadataConfig {
    /// Records per-stage performance measurements when metadata is generated.
    pub track_performance: bool,
    /// Records optimization decisions when metadata is generated.
    pub track_optimization: bool,
    /// Supplies an optional channel diameter for generated metadata.
    pub channel_diameter_m: Option<Length<f64>>,
}

impl MetadataConfig {
    /// Sets the channel diameter used by generated metadata.
    #[must_use]
    pub const fn with_channel_diameter_m(mut self, channel_diameter_m: Length<f64>) -> Self {
        self.channel_diameter_m = Some(channel_diameter_m);
        self
    }
}

/// Generates a network blueprint from geometric dimensions and split stages.
#[must_use]
pub fn create_geometry(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
) -> NetworkBlueprint {
    let total_branches = splits
        .iter()
        .map(SplitType::branch_count)
        .product::<usize>()
        .max(1);
    GeometryGenerator::new(box_dims, *config, *channel_type_config, total_branches).generate(splits)
}

/// Generates a network blueprint through the blueprint-named entry point.
#[must_use]
pub fn create_blueprint_geometry(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
) -> NetworkBlueprint {
    create_geometry(box_dims, splits, config, channel_type_config)
}

/// Generates a network blueprint and attaches configured generation metadata.
#[must_use]
pub fn create_geometry_with_metadata(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
    metadata_config: &MetadataConfig,
) -> NetworkBlueprint {
    let total_branches = splits
        .iter()
        .map(SplitType::branch_count)
        .product::<usize>()
        .max(1);
    GeometryGenerator::new_with_metadata(
        box_dims,
        *config,
        *channel_type_config,
        total_branches,
        metadata_config.clone(),
    )
    .generate(splits)
}

/// Generates a metadata-bearing network blueprint through the blueprint-named entry point.
#[must_use]
pub fn create_blueprint_geometry_with_metadata(
    box_dims: (f64, f64),
    splits: &[SplitType],
    config: &GeometryConfig,
    channel_type_config: &ChannelTypeConfig,
    metadata_config: &MetadataConfig,
) -> NetworkBlueprint {
    create_geometry_with_metadata(
        box_dims,
        splits,
        config,
        channel_type_config,
        metadata_config,
    )
}

/// Configures geometry generation and optional blueprint metadata before building.
pub struct GeometryGeneratorBuilder<'s> {
    box_dims: (f64, f64),
    splits: &'s [SplitType],
    config: GeometryConfig,
    channel_type_config: ChannelTypeConfig,
    metadata_config: Option<MetadataConfig>,
    render_hints: Option<BlueprintRenderHints>,
    topology: Option<BlueprintTopologySpec>,
    lineage: Option<TopologyLineageMetadata>,
    blueprint_metadata: Option<MetadataContainer>,
}

impl<'s> GeometryGeneratorBuilder<'s> {
    /// Creates a builder borrowing the split sequence and copying generation configuration.
    #[must_use]
    pub fn new(
        box_dims: (f64, f64),
        splits: &'s [SplitType],
        config: &GeometryConfig,
        channel_type_config: &ChannelTypeConfig,
    ) -> Self {
        Self {
            box_dims,
            splits,
            config: *config,
            channel_type_config: *channel_type_config,
            metadata_config: None,
            render_hints: None,
            topology: None,
            lineage: None,
            blueprint_metadata: None,
        }
    }

    /// Adds performance and optimization metadata configuration.
    #[must_use]
    pub fn with_metadata_config(mut self, cfg: MetadataConfig) -> Self {
        self.metadata_config = Some(cfg);
        self
    }

    /// Adds rendering hints to the generated blueprint.
    #[must_use]
    pub fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    /// Adds the topology specification used to describe the generated blueprint.
    #[must_use]
    pub fn with_topology_spec(mut self, spec: BlueprintTopologySpec) -> Self {
        self.topology = Some(spec);
        self
    }

    /// Adds lineage metadata describing the source topology.
    #[must_use]
    pub fn with_lineage(mut self, lineage: TopologyLineageMetadata) -> Self {
        self.lineage = Some(lineage);
        self
    }

    /// Adds the metadata container attached to the generated blueprint.
    #[must_use]
    pub fn with_blueprint_metadata(mut self, metadata: MetadataContainer) -> Self {
        self.blueprint_metadata = Some(metadata);
        self
    }

    /// Builds the configured network blueprint.
    #[must_use]
    pub fn build(self) -> NetworkBlueprint {
        let total_branches = self
            .splits
            .iter()
            .map(SplitType::branch_count)
            .product::<usize>()
            .max(1);

        let mut gen = if let Some(metadata_config) = self.metadata_config {
            GeometryGenerator::new_with_metadata(
                self.box_dims,
                self.config,
                self.channel_type_config,
                total_branches,
                metadata_config,
            )
        } else {
            GeometryGenerator::new(
                self.box_dims,
                self.config,
                self.channel_type_config,
                total_branches,
            )
        };

        if let Some(hints) = self.render_hints {
            gen = gen.with_render_hints(hints);
        }
        if let Some(spec) = self.topology {
            gen = gen.with_topology_spec(spec);
        }
        if let Some(lineage) = self.lineage {
            gen = gen.with_lineage(lineage);
        }
        if let Some(metadata) = self.blueprint_metadata {
            gen = gen.with_blueprint_metadata(metadata);
        }

        gen.generate(self.splits)
    }
}
