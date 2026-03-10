mod composite;
mod n_furcation;
mod serpentine;
mod venturi;

use crate::domain::model::NetworkBlueprint;

pub use composite::{
    asymmetric_bifurcation_serpentine_rect, asymmetric_trifurcation_venturi_rect,
    bifurcation_serpentine_rect, bifurcation_trifurcation_venturi_rect, bifurcation_venturi_rect,
    cascade_center_trifurcation_rect, cascade_tri_bi_tri_selective_rect, cell_separation_rect,
    constriction_expansion_array_rect, double_bifurcation_serpentine_rect,
    double_bifurcation_venturi_rect, double_trifurcation_cif_venturi_rect,
    double_trifurcation_venturi_rect, incremental_filtration_tri_bi_rect,
    incremental_filtration_tri_bi_rect_staged, incremental_filtration_tri_bi_rect_staged_remerge,
    n_furcation_serpentine_rect, n_furcation_venturi_rect, parallel_microchannel_array_rect,
    primitive_selective_split_tree_rect, quad_trifurcation_venturi_rect,
    serial_double_venturi_rect, spiral_channel_rect,
    trifurcation_bifurcation_bifurcation_venturi_rect, trifurcation_bifurcation_venturi_rect,
    trifurcation_serpentine_rect, trifurcation_venturi_rect, triple_bifurcation_venturi_rect,
    triple_trifurcation_venturi_rect, venturi_serpentine_rect, CenterSerpentineSpec,
};
pub use n_furcation::{
    bifurcation_rect, n_furcation_rect, symmetric_bifurcation, symmetric_n_furcation,
    symmetric_trifurcation, trifurcation_rect,
};
pub use serpentine::{serpentine_chain, serpentine_rect, serpentine_venturi_rect};
pub use venturi::{venturi_chain, venturi_rect};

pub(crate) fn finalize_preset_blueprint(mut blueprint: NetworkBlueprint) -> NetworkBlueprint {
    crate::visualizations::schematic::materialize_blueprint_layout(&mut blueprint);
    blueprint
}
