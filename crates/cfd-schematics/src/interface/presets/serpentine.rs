#![allow(deprecated)] // NetworkBlueprint::new() is used intentionally; nodes use NodeSpec::new_at() or topological layout is acceptable for linear chains.

use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::topology::presets::{serpentine_bend_venturi_series_spec, serpentine_series_spec};
use crate::BlueprintTopologyFactory;



#[must_use]
pub fn serpentine_chain(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    diameter_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    for i in 1..segments {
        bp.add_node(NodeSpec::new(format!("turn_{i}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    for i in 0..segments {
        let from = if i == 0 {
            "inlet".to_string()
        } else {
            format!("turn_{i}")
        };

        let to = if i + 1 == segments {
            "outlet".to_string()
        } else {
            format!("turn_{}", i + 1)
        };

        let mut spec = ChannelSpec::new_pipe(
            format!("segment_{}", i + 1),
            from,
            to,
            segment_length_m,
            diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        );
        spec.channel_shape = ChannelShape::Serpentine {
            segments,
            bend_radius_m: diameter_m * 0.5,
        };
        bp.add_channel(spec);
    }

    bp
}

/// Rectangular-channel serpentine network (photolithographic millifluidic chips).
///
/// All segments share the same cross-section `width_m × height_m`.
/// Resistance per segment uses the Shah (1978) Poiseuille-number correction for
/// rectangular ducts.  Segment channels are named `segment_1`, `segment_2`, …
#[must_use]
pub fn serpentine_rect(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    let topology = serpentine_series_spec(&name_str, segments, segment_length_m, width_m, height_m);
    BlueprintTopologyFactory::build(&topology).expect("valid serpentine topology spec")
}

/// Rectangular serpentine with venturi constrictions at each U-turn bend.
///
/// # Physics — Venturi Placement at Dean-Flow Apex
///
/// In a curved rectangular channel the centripetal acceleration is:
///
/// ```text
/// a_c = u² / R_bend
/// ```
///
/// where `u` is the mean axial velocity and `R_bend` is the bend radius
/// of curvature.  This centripetal acceleration drives Dean secondary
/// vortices whose intensity is characterised by the Dean number:
///
/// ```text
/// De = Re · √(D_h / (2 · R_bend))
/// ```
///
/// At the apex of each U-turn (minimum `R_bend`), both `a_c` and `De`
/// are maximal.  Placing a venturi throat here exploits:
///
/// 1. **Dean-enhanced cell margination** — the secondary vortices
///    pre-focus cells toward the channel centreline before entering the
///    throat, increasing the fraction of circulating tumour cells (CTCs)
///    that traverse the high-shear cavitation zone.
///
/// 2. **Maximum wall shear** — centripetal acceleration adds to the
///    axial velocity gradient at the outer wall, producing the highest
///    wall-shear-rate site in the entire serpentine.
///
/// # Topology
///
/// ```text
/// inlet → seg_1 → [approach_1 → throat_1 → recovery_1] → seg_2
///       → [approach_2 → throat_2 → recovery_2] → seg_3
///       → …
///       → [approach_{n-1} → throat_{n-1} → recovery_{n-1}] → seg_n → outlet
/// ```
///
/// Each straight segment has `ChannelShape::Serpentine`; each throat is
/// tagged `TherapyZone::CancerTarget` with full `VenturiGeometryMetadata`.
///
/// # Parameters
///
/// - `segments` — number of straight segments (≥ 2; venturis at n−1 bends)
/// - `segment_length_m` — length of each straight segment
/// - `width_m` — channel width (approach / recovery / straight)
/// - `throat_width_m` — constricted width at each venturi throat
/// - `height_m` — channel height (constant — planar chip fabrication)
/// - `throat_length_m` — axial length of the venturi throat section
/// - `bend_radius_m` — serpentine U-turn radius of curvature
#[must_use]
pub fn serpentine_venturi_rect(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    bend_radius_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    let topology = serpentine_bend_venturi_series_spec(
        &name_str,
        segments,
        segment_length_m,
        width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        bend_radius_m,
    );
    BlueprintTopologyFactory::build(&topology).expect("valid venturi topology spec")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::NodeKind;

    #[test]
    fn serpentine_chain_has_correct_segment_count() {
        let bp = serpentine_chain("test", 4, 0.01, 0.001);
        let segs = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segs, 4);
    }

    #[test]
    fn serpentine_rect_has_correct_endpoints() {
        let bp = serpentine_rect("test", 3, 0.01, 0.002, 0.001);
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
        assert_eq!(inlets, 1);
        assert_eq!(outlets, 1);
    }

    #[test]
    fn serpentine_rect_embeds_canonical_topology_metadata() {
        let bp = serpentine_rect("meta-serp", 4, 0.01, 0.002, 0.001);
        let topology = bp
            .topology_spec()
            .expect("serpentine_rect should attach topology metadata");

        assert_eq!(topology.design_name, "meta-serp");
        assert_eq!(topology.treatment_channel_ids().len(), 4);
        assert!(topology.has_serpentine());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }

    #[test]
    fn serpentine_venturi_rect_has_venturi_at_each_bend() {
        let bp = serpentine_venturi_rect(
            "test_sv", 5,     // 5 segments → 4 U-turns → 4 venturis
            0.01,  // 10 mm segments
            0.002, // 2 mm width
            0.001, // 1 mm throat
            0.001, // 1 mm height
            0.003, // 3 mm throat length
            0.001, // 1 mm bend radius
        );

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
        assert_eq!(inlets, 1);
        assert_eq!(outlets, 1);

        let venturi_channels = bp.venturi_channels();
        assert_eq!(
            venturi_channels.len(),
            4,
            "expected 4 venturi throats (one per bend)"
        );

        let segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segments, 5, "expected 5 straight segments");
    }

    #[test]
    fn serpentine_venturi_rect_minimum_segments() {
        let bp = serpentine_venturi_rect("test_min", 1, 0.01, 0.002, 0.001, 0.001, 0.003, 0.001);
        // segments clamped to 2 → 1 venturi
        let venturi_channels = bp.venturi_channels();
        assert_eq!(venturi_channels.len(), 1);

        let segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segments, 2);
    }

    #[test]
    fn serpentine_venturi_rect_embeds_series_topology_metadata() {
        let bp =
            serpentine_venturi_rect("meta-serp-vent", 4, 0.01, 0.002, 0.001, 0.001, 0.003, 0.001);
        let topology = bp
            .topology_spec()
            .expect("serpentine_venturi_rect should attach topology metadata");

        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 3);
        assert_eq!(topology.treatment_channel_ids().len(), 3);
        assert!(topology.has_serpentine());
    }

    #[test]
    fn serpentine_venturi_rect_throat_geometry_is_correct() {
        let bp = serpentine_venturi_rect("test_geom", 3, 0.01, 0.002, 0.001, 0.001, 0.003, 0.001);

        for ch in bp.venturi_channels() {
            let vg = ch
                .venturi_geometry
                .as_ref()
                .expect("venturi throat must have geometry");
            assert!(
                (vg.throat_width_m - 0.001).abs() < 1e-12,
                "throat width mismatch"
            );
            assert!(
                (vg.inlet_width_m - 0.002).abs() < 1e-12,
                "inlet width mismatch"
            );
            assert!(
                (vg.outlet_width_m - 0.002).abs() < 1e-12,
                "outlet width mismatch"
            );
            assert!(vg.convergent_half_angle_deg > 0.0 && vg.convergent_half_angle_deg < 90.0);
        }
    }
}
