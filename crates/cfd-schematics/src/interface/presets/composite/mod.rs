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
mod series;
mod specialized;
mod trifurcation;

pub use bifurcation::{bifurcation_serpentine_rect, bifurcation_venturi_rect};
pub use mixed_tree::{
    bifurcation_trifurcation_venturi_rect, double_trifurcation_venturi_rect,
    trifurcation_bifurcation_bifurcation_venturi_rect, trifurcation_bifurcation_venturi_rect,
};
pub use multi_level::{
    double_bifurcation_venturi_rect, quad_trifurcation_venturi_rect,
    triple_bifurcation_venturi_rect, triple_trifurcation_venturi_rect,
};
pub use series::{serial_double_venturi_rect, venturi_serpentine_rect};
pub use specialized::{
    asymmetric_bifurcation_serpentine_rect, asymmetric_trifurcation_venturi_rect,
    cascade_center_trifurcation_rect, cascade_tri_bi_tri_selective_rect, cell_separation_rect,
    constriction_expansion_array_rect, incremental_filtration_tri_bi_rect,
    incremental_filtration_tri_bi_rect_staged, incremental_filtration_tri_bi_rect_staged_remerge,
    parallel_microchannel_array_rect, spiral_channel_rect,
};
pub use trifurcation::{trifurcation_serpentine_rect, trifurcation_venturi_rect};

pub(crate) const BLOOD_MU: f64 = 3.5e-3; // Pa·s (whole blood, high-shear Newtonian approx.)

/// Shah-London hydraulic resistance for a rectangular duct [Pa·s/m³].
///
/// # Theorem (Shah & London 1978)
///
/// For fully-developed laminar flow through a rectangular duct with aspect
/// ratio α = min(w,h)/max(w,h), the product fRe (Poiseuille number) is:
///
/// ```text
/// Po = 96 · (1 − 1.3553α + 1.9467α² − 1.7012α³ + 0.9564α⁴ − 0.2537α⁵)
/// ```
///
/// The hydraulic resistance is then `R = Po·μ·L / (D_h² · A)`.
pub(crate) fn shah_london(w: f64, h: f64, l: f64, mu: f64) -> f64 {
    let alpha = h.min(w) / h.max(w);
    let po = 96.0
        * (1.0 - 1.3553 * alpha + 1.9467 * alpha.powi(2) - 1.7012 * alpha.powi(3)
            + 0.9564 * alpha.powi(4)
            - 0.2537 * alpha.powi(5));
    let d_h = 2.0 * w * h / (w + h);
    let area = w * h;
    po * mu * l / (d_h * d_h * area)
}

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
    fn bifurcation_venturi_has_single_inlet_outlet() {
        let bp = bifurcation_venturi_rect("t", 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn trifurcation_venturi_has_single_inlet_outlet() {
        let bp = trifurcation_venturi_rect("t", 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn bifurcation_serpentine_has_single_inlet_outlet() {
        let bp = bifurcation_serpentine_rect("t", 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn trifurcation_serpentine_has_single_inlet_outlet() {
        let bp = trifurcation_serpentine_rect("t", 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
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
