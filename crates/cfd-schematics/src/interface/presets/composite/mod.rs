//! Composite closed-loop preset networks for millifluidic SDT device designs.
//!
//! Every function in this module produces a [`NetworkBlueprint`] with exactly
//! **one `NodeKind::Inlet`** and **one `NodeKind::Outlet`**.  Every split
//! junction is paired with a corresponding converging junction so the full
//! 127.76 × 85.47 mm 96-well-plate cuboid is always traversed completely.
//!
//! # Channel-name conventions (relied on by downstream consumers)
//! - `"inlet_section"` — first inlet / trunk approach channel
//! - `"throat_section"` — primary venturi throat (may also be `"throat_section_N"`)
//! - `"segment_1"`, `"segment_2"`, … — serpentine straight segments
//! - `"trunk_out"` — final trunk channel before the outlet

mod bifurcation;
mod mixed_tree;
mod multi_level;
mod n_furcation;
mod series;
mod specialized;
mod trifurcation;

pub use bifurcation::{
    bifurcation_serpentine_rect, bifurcation_venturi_rect, double_bifurcation_serpentine_rect,
};
pub use mixed_tree::{
    bifurcation_trifurcation_venturi_rect, double_trifurcation_venturi_rect,
    trifurcation_bifurcation_bifurcation_venturi_rect, trifurcation_bifurcation_venturi_rect,
};
pub use multi_level::{
    double_bifurcation_venturi_rect, quad_trifurcation_venturi_rect,
    triple_bifurcation_venturi_rect, triple_trifurcation_venturi_rect,
};
pub use n_furcation::{n_furcation_serpentine_rect, n_furcation_venturi_rect};
pub use series::{serial_double_venturi_rect, venturi_serpentine_rect};
pub use specialized::{
    asymmetric_bifurcation_serpentine_rect, asymmetric_trifurcation_venturi_rect,
    cascade_center_trifurcation_rect, cascade_tri_bi_tri_selective_rect, cell_separation_rect,
    constriction_expansion_array_rect, double_trifurcation_cif_venturi_rect,
    incremental_filtration_tri_bi_rect, incremental_filtration_tri_bi_rect_staged,
    incremental_filtration_tri_bi_rect_staged_remerge, parallel_microchannel_array_rect,
    primitive_selective_split_tree_rect, spiral_channel_rect, CenterSerpentineSpec,
};
pub use trifurcation::{trifurcation_serpentine_rect, trifurcation_venturi_rect};

pub(crate) use super::finalize_preset_blueprint;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::{NetworkBlueprint, NodeKind};

    fn count_inlets_outlets(bp: &NetworkBlueprint) -> (usize, usize) {
        let inlets = bp
            .nodes
            .iter()
            .filter(|n| n.kind == NodeKind::Inlet)
            .count();
        let outlets = bp
            .nodes
            .iter()
            .filter(|n| n.kind == NodeKind::Outlet)
            .count();
        (inlets, outlets)
    }

    #[test]
    fn venturi_serpentine_has_single_inlet_outlet() {
        let bp = venturi_serpentine_rect("t", 2e-3, 0.5e-3, 0.5e-3, 1e-3, 6, 7.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
    }

    #[test]
    fn n_furcation_venturi_has_single_inlet_outlet() {
        let bp = n_furcation_venturi_rect("t", 2, 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn n_furcation_serpentine_has_single_inlet_outlet() {
        let bp = n_furcation_serpentine_rect("t", 3, 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn cell_separation_has_single_inlet_outlet() {
        let bp = cell_separation_rect("t", 22.5e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3, 22.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }
}
