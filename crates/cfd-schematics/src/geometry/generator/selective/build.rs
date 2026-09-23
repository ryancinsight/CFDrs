use super::builder::SelectiveTreeBuilder;
use super::request::{SelectiveTreeRequest, SelectiveTreeTopology};
use crate::domain::model::NetworkBlueprint;

/// Generates a selective-tree blueprint from its physical request and topology.
pub fn create_selective_tree_geometry(request: &SelectiveTreeRequest) -> NetworkBlueprint {
    let geometry = request.geometry();
    let mut builder = SelectiveTreeBuilder::new(request.name.clone(), geometry.box_dims_mm);
    match &request.topology {
        SelectiveTreeTopology::CascadeCenterTrifurcation {
            n_levels,
            center_frac,
            venturi_treatment_enabled,
            center_serpentine,
        } => builder.build_cct(
            *n_levels,
            *center_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            &geometry,
        ),
        SelectiveTreeTopology::IncrementalFiltrationTriBi {
            n_pretri,
            pretri_center_frac,
            terminal_tri_center_frac,
            bi_treat_frac,
            venturi_treatment_enabled,
            center_serpentine,
            outlet_tail_length_m,
        } => builder.build_cif(
            *n_pretri,
            *pretri_center_frac,
            *terminal_tri_center_frac,
            *bi_treat_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            outlet_tail_length_m.into_base(),
            &geometry,
        ),
        SelectiveTreeTopology::TriBiTriSelective {
            first_center_frac,
            bi_treat_frac,
            second_center_frac,
        } => builder.build_tbt(
            *first_center_frac,
            *bi_treat_frac,
            *second_center_frac,
            &geometry,
        ),
        SelectiveTreeTopology::DoubleTrifurcationCif {
            split1_center_frac,
            split2_center_frac,
            center_throat_count,
            inter_throat_spacing_m,
        } => builder.build_dtcv(
            *split1_center_frac,
            *split2_center_frac,
            *center_throat_count,
            inter_throat_spacing_m.into_base(),
            &geometry,
        ),
    }
    builder.finish()
}

