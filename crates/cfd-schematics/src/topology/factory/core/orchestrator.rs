use crate::domain::model::NetworkBlueprint;
use crate::error::Result;
use crate::topology::model::BlueprintTopologySpec;

/// Core interface for turning declarative [`BlueprintTopologySpec`]s into
/// [`NetworkBlueprint`] graphs.
///
/// ## SSOT Architecture
///
/// This factory is a **thin facade** that delegates all geometry generation
/// to the canonical
/// [`GeometryGeneratorBuilder`](crate::geometry::generator::GeometryGeneratorBuilder)
/// pipeline. No ad-hoc
/// node/channel construction is performed here — that logic lives exclusively
/// in the private `GeometryGenerator` implementation.
pub struct BlueprintTopologyFactory;

impl BlueprintTopologyFactory {
    /// Entrypoint: builds a fully detailed blueprint graph from a declarative
    /// topology spec by delegating to the canonical `create_geometry` pipeline.
    ///
    /// # Pipeline
    ///
    /// 1. Validate the spec (`validation::validate_spec`)
    /// 2. Convert `BlueprintTopologySpec` → `SplitType[]` + `GeometryConfig`
    /// 3. Delegate to `GeometryGeneratorBuilder` (canonical pipeline)
    /// 4. Apply venturi placements post-hoc
    /// 5. Attach topology + lineage metadata
    ///
    /// # Errors
    ///
    /// Returns a descriptive error string if the spec violates any geometric
    /// or structural constraint.
    pub fn build(spec: &BlueprintTopologySpec) -> Result<NetworkBlueprint> {
        super::super::validation::validate_spec(spec)?;

        let lineage = Self::lineage_for_spec(spec);

        let mut blueprint = if spec.has_series_path() && spec.split_stages.is_empty() {
            Self::build_series_path(spec, lineage)
        } else if spec.has_parallel_paths() && spec.split_stages.is_empty() {
            Self::build_parallel_path(spec, lineage)
        } else {
            Self::build_split_tree(spec, lineage)?
        };
        let resolved_spec = Self::resolve_materialized_venturi_targets(&blueprint, spec);
        blueprint.topology = Some(resolved_spec.clone());

        // Post-process: apply venturi placements
        if resolved_spec.has_venturi()
            && !Self::has_materialized_venturi_geometry(&blueprint, &resolved_spec)
        {
            super::super::modifiers::venturi::apply_venturi_placements(&mut blueprint, &resolved_spec)?;
        }

        Ok(blueprint)
    }

    // build_impl.rs: build_series_path, build_parallel_path, build_split_tree,
    // reconcile_channel_ids, venturi helpers, leading_merge_side_treatment_channels

    // mutation_impl.rs: validate_spec, mutate

    // spec_analysis_impl.rs: estimate_dean_site, lineage_for_spec, spec queries
}
