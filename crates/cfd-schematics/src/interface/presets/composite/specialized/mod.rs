//! Specialized composite presets: cell separation, asymmetric, constriction,
//! spiral, parallel microchannel array, and selective tree topologies.

mod basic;
mod filtration;
mod parallel_lane;
mod trifurcation;

pub use basic::{
    asymmetric_bifurcation_serpentine_rect, cell_separation_rect,
    constriction_expansion_array_rect, parallel_microchannel_array_rect,
    primitive_selective_split_tree_rect, spiral_channel_rect,
};
pub use filtration::{
    cascade_center_trifurcation_rect, incremental_filtration_tri_bi_rect,
    incremental_filtration_tri_bi_rect_staged, incremental_filtration_tri_bi_rect_staged_remerge,
};
pub use parallel_lane::CenterSerpentineSpec;
pub use trifurcation::{
    asymmetric_trifurcation_venturi_rect, cascade_tri_bi_tri_selective_rect,
    double_trifurcation_cif_venturi_rect,
};

// Sub-modules provide their own imports

#[cfg(test)]
mod tests;
