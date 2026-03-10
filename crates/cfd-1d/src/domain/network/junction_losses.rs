//! Family-aware minor-loss augmentation for blueprint junctions.
//!
//! The 1D network solve remains edge-based, but junctions with explicit
//! geometry metadata receive family-specific quadratic loss terms. This keeps
//! the solver fast while making bifurcation/trifurcation/merge behavior more
//! faithful than a generic T-junction heuristic.

use super::wrapper::Network;
use cfd_core::physics::fluid::FluidTrait;
use cfd_schematics::domain::model::NetworkBlueprint;
use cfd_schematics::geometry::metadata::{
    JunctionFamily, JunctionGeometryMetadata, MetadataContainer,
};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct JunctionLossAuditStats {
    pub explicit_metadata_nodes: usize,
    pub fallback_nodes: usize,
}

#[derive(Debug, Clone, Copy)]
struct JunctionLossProfile {
    branch_k: f64,
    run_k: f64,
    area_ratio_exp: f64,
    branch_angle_ref_deg: f64,
    run_angle_ref_deg: f64,
}

pub(crate) fn apply_blueprint_junction_losses<T, F>(
    network: &mut Network<T, F>,
    blueprint: &NetworkBlueprint,
) -> JunctionLossAuditStats
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    use crate::domain::network::NodeType;

    let mut stats = JunctionLossAuditStats::default();
    // Derive fluid density from the network's fluid model at reference conditions
    // (body temperature 310.15 K, atmospheric pressure 101 325 Pa).
    // Falls back to 1060 kg/m³ (whole-blood reference) if the fluid query fails.
    let rho_ref_t = T::from_f64(310.15).expect("Mathematical constant conversion compromised");
    let p_ref_t = T::from_f64(101_325.0).expect("Mathematical constant conversion compromised");
    let rho_blood: T = network
        .fluid
        .properties_at(rho_ref_t, p_ref_t)
        .map(|state| state.density)
        .unwrap_or_else(|_| {
            T::from_f64(1060.0).expect("Mathematical constant conversion compromised")
        });
    let two: T = T::one() + T::one();

    let blueprint_node_meta: HashMap<
        &str,
        (
            Option<&JunctionGeometryMetadata>,
            Option<&MetadataContainer>,
        ),
    > = blueprint
        .nodes
        .iter()
        .map(|node| {
            (
                node.id.as_str(),
                (node.junction_geometry.as_ref(), node.metadata.as_ref()),
            )
        })
        .collect();

    let junction_info: Vec<(petgraph::graph::NodeIndex, usize, usize, usize)> = network
        .graph
        .node_indices()
        .filter_map(|idx| {
            let node = network.graph.node_weight(idx)?;
            if node.node_type != NodeType::Junction {
                return None;
            }
            let in_deg = network
                .graph
                .edges_directed(idx, petgraph::Direction::Incoming)
                .count();
            let out_deg = network
                .graph
                .edges_directed(idx, petgraph::Direction::Outgoing)
                .count();
            let total = in_deg + out_deg;
            (total >= 3).then_some((idx, in_deg, out_deg, total))
        })
        .collect();

    for (node_idx, in_deg, out_deg, total) in junction_info {
        let outgoing: Vec<petgraph::graph::EdgeIndex> = network
            .graph
            .edges_directed(node_idx, petgraph::Direction::Outgoing)
            .map(|e| e.id())
            .collect();
        let incoming: Vec<petgraph::graph::EdgeIndex> = network
            .graph
            .edges_directed(node_idx, petgraph::Direction::Incoming)
            .map(|e| e.id())
            .collect();

        let junction_meta = network
            .graph
            .node_weight(node_idx)
            .map(|node| node.id.as_str())
            .and_then(|node_name| blueprint_node_meta.get(node_name).copied())
            .and_then(|(geometry, metadata)| {
                geometry
                    .or_else(|| metadata.and_then(|meta| meta.get::<JunctionGeometryMetadata>()))
            });

        if junction_meta.is_some() {
            stats.explicit_metadata_nodes += 1;
        } else {
            stats.fallback_nodes += 1;
        }

        let family = junction_meta.map(|meta| meta.junction_family);
        let profile = junction_loss_profile(family, in_deg, out_deg, total);
        let (run_edges, branch_edges) =
            classify_run_and_branch_edges(family, in_deg, out_deg, &incoming, &outgoing);

        let branch_angle = junction_meta
            .and_then(|meta| representative_abs_angle_deg(&meta.branch_angles_deg))
            .or_else(|| {
                junction_meta.and_then(|meta| representative_abs_angle_deg(&meta.merge_angles_deg))
            });
        let run_angle = junction_meta
            .and_then(|meta| representative_abs_angle_deg(&meta.merge_angles_deg))
            .or_else(|| {
                junction_meta.and_then(|meta| representative_abs_angle_deg(&meta.branch_angles_deg))
            });

        let k_branch = scale_k_by_angle(
            profile.branch_k,
            branch_angle,
            profile.branch_angle_ref_deg,
            false,
        );
        let k_run = scale_k_by_angle(profile.run_k, run_angle, profile.run_angle_ref_deg, true);

        let k_branch_t =
            T::from_f64(k_branch).expect("Mathematical constant conversion compromised");
        let k_run_t = T::from_f64(k_run).expect("Mathematical constant conversion compromised");
        let exp_t = T::from_f64(profile.area_ratio_exp)
            .expect("Mathematical constant conversion compromised");

        let run_area = run_edges
            .first()
            .and_then(|eidx| network.graph.edge_weight(*eidx))
            .map(|edge| edge.area);
        let branch_area_sum = branch_edges.iter().fold(T::zero(), |acc, eidx| {
            acc + network
                .graph
                .edge_weight(*eidx)
                .map_or(T::zero(), |edge| edge.area)
        });

        for eidx in &branch_edges {
            if let Some(edge) = network.graph.edge_weight_mut(*eidx) {
                let area_sq = edge.area * edge.area;
                if area_sq <= T::zero() {
                    continue;
                }
                let area_ratio_scale = run_area
                    .filter(|area| *area > T::zero())
                    .map_or(T::one(), |run_area| (edge.area / run_area).powf(exp_t));
                let k_correction = (k_branch_t * area_ratio_scale) * rho_blood / (two * area_sq);
                edge.quad_coeff += k_correction;
            }
        }

        for eidx in &run_edges {
            if let Some(edge) = network.graph.edge_weight_mut(*eidx) {
                let area_sq = edge.area * edge.area;
                if area_sq <= T::zero() {
                    continue;
                }
                let area_ratio_scale = if branch_area_sum > T::zero() && edge.area > T::zero() {
                    (branch_area_sum / edge.area).powf(
                        T::from_f64(0.25).expect("Mathematical constant conversion compromised"),
                    )
                } else {
                    T::one()
                };
                let k_correction = (k_run_t * area_ratio_scale) * rho_blood / (two * area_sq);
                edge.quad_coeff += k_correction;
            }
        }
    }

    stats
}

fn classify_run_and_branch_edges(
    family: Option<JunctionFamily>,
    in_deg: usize,
    out_deg: usize,
    incoming: &[petgraph::graph::EdgeIndex],
    outgoing: &[petgraph::graph::EdgeIndex],
) -> (
    Vec<petgraph::graph::EdgeIndex>,
    Vec<petgraph::graph::EdgeIndex>,
) {
    match family {
        Some(JunctionFamily::Cross) => {
            (Vec::new(), [incoming.to_vec(), outgoing.to_vec()].concat())
        }
        Some(JunctionFamily::Merge) => (outgoing.to_vec(), incoming.to_vec()),
        Some(JunctionFamily::Bifurcation | JunctionFamily::Trifurcation | JunctionFamily::Tee) => {
            if in_deg == 1 && out_deg >= 2 {
                (incoming.to_vec(), outgoing.to_vec())
            } else if in_deg >= 2 && out_deg == 1 {
                (outgoing.to_vec(), incoming.to_vec())
            } else {
                (Vec::new(), [incoming.to_vec(), outgoing.to_vec()].concat())
            }
        }
        None => {
            if in_deg == 1 && out_deg >= 2 {
                (incoming.to_vec(), outgoing.to_vec())
            } else if in_deg >= 2 && out_deg == 1 {
                (outgoing.to_vec(), incoming.to_vec())
            } else {
                (Vec::new(), [incoming.to_vec(), outgoing.to_vec()].concat())
            }
        }
    }
}

fn junction_loss_profile(
    family: Option<JunctionFamily>,
    in_deg: usize,
    out_deg: usize,
    total: usize,
) -> JunctionLossProfile {
    match family {
        Some(JunctionFamily::Bifurcation) if in_deg == 1 && out_deg >= 2 => JunctionLossProfile {
            branch_k: 0.90,
            run_k: 0.03,
            area_ratio_exp: 0.45,
            branch_angle_ref_deg: 35.0,
            run_angle_ref_deg: 5.0,
        },
        Some(JunctionFamily::Trifurcation) if in_deg == 1 && out_deg >= 3 => JunctionLossProfile {
            branch_k: 1.05,
            run_k: 0.04,
            area_ratio_exp: 0.40,
            branch_angle_ref_deg: 25.0,
            run_angle_ref_deg: 5.0,
        },
        Some(JunctionFamily::Merge) => JunctionLossProfile {
            branch_k: 1.35,
            run_k: 0.06,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 30.0,
            run_angle_ref_deg: 10.0,
        },
        Some(JunctionFamily::Tee) if in_deg == 1 && out_deg >= 2 => JunctionLossProfile {
            branch_k: 1.00,
            run_k: 0.05,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 90.0,
            run_angle_ref_deg: 0.0,
        },
        Some(JunctionFamily::Tee) => JunctionLossProfile {
            branch_k: 1.50,
            run_k: 0.04,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 45.0,
            run_angle_ref_deg: 0.0,
        },
        Some(JunctionFamily::Cross) if total == 4 && in_deg == 2 && out_deg == 2 => {
            JunctionLossProfile {
                branch_k: 1.20,
                run_k: 0.05,
                area_ratio_exp: 0.50,
                branch_angle_ref_deg: 90.0,
                run_angle_ref_deg: 90.0,
            }
        }
        Some(JunctionFamily::Cross) => JunctionLossProfile {
            branch_k: 1.30,
            run_k: 0.05,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 90.0,
            run_angle_ref_deg: 90.0,
        },
        Some(JunctionFamily::Bifurcation | JunctionFamily::Trifurcation)
            if in_deg >= 2 && out_deg == 1 =>
        {
            JunctionLossProfile {
                branch_k: 1.30,
                run_k: 0.06,
                area_ratio_exp: 0.45,
                branch_angle_ref_deg: 25.0,
                run_angle_ref_deg: 10.0,
            }
        }
        _ if total == 4 && in_deg == 2 && out_deg == 2 => JunctionLossProfile {
            branch_k: 1.30,
            run_k: 0.05,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 90.0,
            run_angle_ref_deg: 90.0,
        },
        _ if in_deg == 1 && out_deg >= 2 => JunctionLossProfile {
            branch_k: 1.00,
            run_k: 0.05,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 90.0,
            run_angle_ref_deg: 0.0,
        },
        _ if in_deg >= 2 && out_deg == 1 => JunctionLossProfile {
            branch_k: 1.50,
            run_k: 0.04,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 45.0,
            run_angle_ref_deg: 0.0,
        },
        _ => JunctionLossProfile {
            branch_k: 1.00,
            run_k: 0.05,
            area_ratio_exp: 0.50,
            branch_angle_ref_deg: 45.0,
            run_angle_ref_deg: 0.0,
        },
    }
}

fn representative_abs_angle_deg(angles: &[f64]) -> Option<f64> {
    angles
        .iter()
        .copied()
        .filter(|angle| angle.is_finite())
        .map(f64::abs)
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

fn scale_k_by_angle(
    base_k: f64,
    angle_deg: Option<f64>,
    reference_deg: f64,
    run_edge: bool,
) -> f64 {
    let Some(angle_deg) = angle_deg else {
        return base_k;
    };
    let reference = reference_deg.max(15.0);
    let deviation = ((angle_deg.abs() - reference).abs() / reference).clamp(0.0, 2.0);
    let scale = if run_edge {
        1.0 + 0.10 * deviation
    } else {
        1.0 + 0.25 * deviation
    };
    (base_k * scale).max(0.0)
}
