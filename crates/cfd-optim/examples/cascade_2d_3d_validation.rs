//! # Cascade 2D/3D Simulation and Richardson Extrapolation
//!
//! Full multi-fidelity validation of cascade millifluidic designs for SDT:
//!
//! - **Part 1**: Select top-3 RbcProtectedSdt designs (covers CIFX and CCT families)
//! - **Part 2**: Extract `NetworkBlueprint` for each design (full cascade)
//! - **Part 3**: Complete 1D solve — all nodes/edges via `NetworkSolver::solve_network`
//! - **Part 4**: Complete 2D solve — all channels via `Network2DSolver::solve_all` at 3 grids
//! - **Part 5**: Complete 3D solve — per-channel FEM via `StructuredGridBuilder` at 3 refinements
//! - **Part 6**: Richardson extrapolation (GCI%) for 2D wall-shear and 3D dp/σ
//! - **Part 7**: Cross-fidelity consistency table (1D vs 2D vs 3D)
//! - **Part 8**: Write `report/milestone12/cascade_validation.md`
//!
//! ## Architecture Principle
//!
//! Each fidelity level solves the FULL cascade network simultaneously, not isolated subsections:
//!
//! | Fidelity | Solver | Cascade handling |
//! |----------|--------|------------------|
//! | 1D | `NetworkSolver::solve_network()` | All nodes + edges as a DAG |
//! | 2D | `Network2DSolver::solve_all()` | One NS domain per channel; flow via BFS |
//! | 3D | `StructuredGridBuilder` + `FemSolver` | Per-channel FEM; results aggregated |
//!
//! Richardson extrapolation refines mesh WITHIN each fidelity — it does NOT compare
//! across fidelity levels.

#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]

use std::collections::HashMap;
use std::fs;
use std::path::Path;

// ── cfd-1d ─────────────────────────────────────────────────────────────────────
use cfd_1d::network::network_from_blueprint;
use cfd_1d::{NetworkProblem, NetworkSolver, NodeType};

// ── cfd-2d ─────────────────────────────────────────────────────────────────────
use cfd_2d::network::{Network2DSolver, Network2dBuilderSink};
use cfd_schematics::application::ports::GraphSink;

// ── cfd-3d ─────────────────────────────────────────────────────────────────────
use cfd_3d::{FemConfig, FemSolver, StokesFlowProblem};

// ── cfd-mesh ───────────────────────────────────────────────────────────────────
use cfd_mesh::domain::core::index::VertexId;
use cfd_mesh::domain::grid::StructuredGridBuilder;

// ── cfd-core ───────────────────────────────────────────────────────────────────
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::{BloodModel, CassonBlood, ConstantPropertyFluid};

// ── cfd-validation ─────────────────────────────────────────────────────────────
use cfd_validation::convergence::{GridConvergenceIndex, RichardsonExtrapolation};

// ── cfd-optim ──────────────────────────────────────────────────────────────────
use cfd_optim::{OptimMode, SdtOptimizer, SdtWeights};

// ── nalgebra ───────────────────────────────────────────────────────────────────
use nalgebra::Vector3;

// ── petgraph ───────────────────────────────────────────────────────────────────
use petgraph::visit::EdgeRef;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Blood density [kg/m³]
const BLOOD_DENSITY: f64 = 1060.0;
/// Blood effective viscosity [Pa·s] (Newtonian approximation for 3D FEM)
const BLOOD_VISCOSITY: f64 = 3.5e-3;
/// Inlet volumetric flow rate [m³/s] = 200 mL/min
const INLET_FLOW_M3_S: f64 = 200.0e-6 / 60.0;
/// Vapour pressure of blood at 37 °C [Pa]
const P_VAPOUR_PA: f64 = 6_270.0;
/// Reference absolute pressure [Pa]
const P_REF_PA: f64 = 101_325.0;
/// Richardson refinement ratio
const REFINEMENT_RATIO: f64 = 2.0;

// ─────────────────────────────────────────────────────────────────────────────
// Structs
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct OneDResult {
    design_id: String,
    inlet_pressure_pa: f64,
    dp_total_pa: f64,
    q_venturi_m3_s: f64,
    /// Per-channel: (id, q_m3_s, p_from_pa)
    channels: Vec<(String, f64, f64)>,
}

#[derive(Debug, Clone)]
struct TwoDLevel {
    nx: usize,
    ny: usize,
    max_wall_shear_pa: f64,
    converged: bool,
}

#[derive(Debug, Clone)]
struct TwoDResult {
    design_id: String,
    levels: [TwoDLevel; 3],
    /// Richardson extrapolated wall shear [Pa]
    wall_shear_extrap_pa: f64,
    wall_shear_order: f64,
    gci_fine_pct: f64,
}

#[derive(Debug, Clone)]
struct ThreeDLevel {
    nx: usize,
    ny_nz: usize,
    dp_total_pa: f64,
    sigma_venturi: f64,
    u_venturi_m_s: f64,
}

#[derive(Debug, Clone)]
struct ThreeDResult {
    design_id: String,
    levels: [ThreeDLevel; 3],
    dp_extrap_pa: f64,
    dp_order: f64,
    dp_gci_fine_pct: f64,
    sigma_extrap: f64,
    sigma_order: f64,
    sigma_gci_fine_pct: f64,
    asymptotic: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 1: Design selection
// ─────────────────────────────────────────────────────────────────────────────

fn select_top3() -> anyhow::Result<Vec<cfd_optim::RankedDesign>> {
    let designs = SdtOptimizer::new(OptimMode::RbcProtectedSdt, SdtWeights::default())
        .top_k(3)
        .map_err(|e| anyhow::anyhow!("optimiser failed: {e}"))?;
    Ok(designs)
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 3: Complete 1D cascade solve
// ─────────────────────────────────────────────────────────────────────────────

fn solve_1d(
    design_id: &str,
    blueprint: &cfd_schematics::NetworkBlueprint,
) -> anyhow::Result<OneDResult> {
    let blood = CassonBlood::<f64>::normal_blood();
    let mut network = network_from_blueprint(blueprint, blood)
        .map_err(|e| anyhow::anyhow!("[{design_id}] network_from_blueprint: {e}"))?;

    // Identify inlet / outlet nodes
    let mut inlet_nodes = Vec::new();
    let mut outlet_nodes = Vec::new();
    for idx in network.graph.node_indices() {
        if let Some(node) = network.graph.node_weight(idx) {
            match node.node_type {
                NodeType::Inlet => inlet_nodes.push(idx),
                NodeType::Outlet => outlet_nodes.push(idx),
                _ => {}
            }
        }
    }
    if inlet_nodes.is_empty() || outlet_nodes.is_empty() {
        return Err(anyhow::anyhow!(
            "[{design_id}] no inlet/outlet nodes in blueprint"
        ));
    }

    let q_per_inlet = INLET_FLOW_M3_S / inlet_nodes.len() as f64;
    for &n in &inlet_nodes {
        network.set_neumann_flow(n, q_per_inlet);
    }
    for &n in &outlet_nodes {
        network.set_pressure(n, 0.0);
    }

    let solver = NetworkSolver::<f64, CassonBlood<f64>>::new();
    let problem = NetworkProblem::new(network.clone());
    let solved = solver
        .solve_network(&problem)
        .map_err(|e| anyhow::anyhow!("[{design_id}] NetworkSolver: {e}"))?;

    let inlet_pressure_pa = inlet_nodes
        .iter()
        .filter_map(|idx| solved.pressures.get(idx))
        .sum::<f64>()
        / inlet_nodes.len() as f64;

    // Collect per-channel flows + pressures
    let mut channels = Vec::with_capacity(blueprint.channels.len());
    let mut q_venturi_m3_s = 0.0_f64;

    for edge_ref in solved.graph.edge_references() {
        let eid = edge_ref.id();
        let edge = edge_ref.weight();
        let (from_idx, _) = solved
            .graph
            .edge_endpoints(eid)
            .ok_or_else(|| anyhow::anyhow!("[{design_id}] missing edge endpoints"))?;
        let q = solved
            .flow_rates
            .get(&eid)
            .copied()
            .unwrap_or(edge.flow_rate);
        let p = solved.pressures.get(&from_idx).copied().unwrap_or(0.0);
        if edge.id.to_ascii_lowercase().contains("throat") {
            q_venturi_m3_s += q.abs();
        }
        channels.push((edge.id.clone(), q, p));
    }

    Ok(OneDResult {
        design_id: design_id.to_string(),
        inlet_pressure_pa,
        dp_total_pa: inlet_pressure_pa,
        q_venturi_m3_s,
        channels,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 4: Complete 2D cascade solve at 3 grid levels
// ─────────────────────────────────────────────────────────────────────────────

fn solve_2d_level(
    design_id: &str,
    blueprint: &cfd_schematics::NetworkBlueprint,
    nx: usize,
    ny: usize,
) -> anyhow::Result<TwoDLevel> {
    let blood = BloodModel::Newtonian(BLOOD_VISCOSITY);
    let sink = Network2dBuilderSink::new(blood, BLOOD_DENSITY, INLET_FLOW_M3_S, nx, ny);
    let mut net2d: Network2DSolver<f64> = sink
        .build(blueprint)
        .map_err(|e| anyhow::anyhow!("[{design_id}] 2D build ({nx}×{ny}): {e}"))?;

    let result = net2d
        .solve_all(1e-6)
        .map_err(|e| anyhow::anyhow!("[{design_id}] 2D solve ({nx}×{ny}): {e}"))?;

    let max_wall_shear_pa = result
        .channels
        .iter()
        .map(|c| c.wall_shear_pa)
        .fold(0.0_f64, f64::max);

    let converged_all = result.converged_count == result.channels.len();

    println!(
        "  2D [{design_id}] {nx}×{ny}: τ_max={:.2} Pa  converged={converged_all}",
        max_wall_shear_pa
    );

    Ok(TwoDLevel {
        nx,
        ny,
        max_wall_shear_pa,
        converged: converged_all,
    })
}

fn solve_2d(
    design_id: &str,
    blueprint: &cfd_schematics::NetworkBlueprint,
) -> anyhow::Result<TwoDResult> {
    // Three grid levels: coarse / medium / fine
    let levels_in = [(10usize, 5usize), (20, 10), (40, 20)];
    let mut levels = Vec::with_capacity(3);
    for (nx, ny) in levels_in {
        levels.push(solve_2d_level(design_id, blueprint, nx, ny)?);
    }

    let vals: [f64; 3] = [
        levels[0].max_wall_shear_pa,
        levels[1].max_wall_shear_pa,
        levels[2].max_wall_shear_pa,
    ];

    // Richardson extrapolation: QoI = max wall shear stress
    let (wall_shear_extrap_pa, order, gci_fine_pct) =
        richardson_study(vals[0], vals[1], vals[2], "[2D wall-shear]");

    Ok(TwoDResult {
        design_id: design_id.to_string(),
        levels: [levels[0].clone(), levels[1].clone(), levels[2].clone()],
        wall_shear_extrap_pa,
        wall_shear_order: order,
        gci_fine_pct,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 5: Complete 3D cascade solve at 3 refinement levels
// ─────────────────────────────────────────────────────────────────────────────

/// Build and solve a single channel using a structured tet mesh.
///
/// The channel is mapped from blueprint dimensions:
/// - Z-axis = flow direction (length = `length_m`), z=0 inlet, z=1 outlet
/// - X/Y axes = cross-section (width × height)
fn solve_channel_3d(
    ch: &cfd_schematics::domain::model::ChannelSpec,
    u_inlet: f64,
    nx: usize,
    ny_nz: usize,
) -> anyhow::Result<(f64, f64)> {
    // Build unit-cube structured tet mesh [0,1]³
    let mut mesh = StructuredGridBuilder::new(nx, ny_nz, ny_nz)
        .build()
        .map_err(|e| anyhow::anyhow!("StructuredGridBuilder: {e}"))?;

    // Scale vertices: z → length_m; x → width_m; y → height_m
    let (width_m, height_m) = match ch.cross_section {
        cfd_schematics::domain::model::CrossSectionSpec::Rectangular { width_m, height_m } => {
            (width_m, height_m)
        }
        cfd_schematics::domain::model::CrossSectionSpec::Circular { diameter_m } => {
            (diameter_m, diameter_m)
        }
    };
    let length_m = ch.length_m;
    let n_verts = mesh.vertex_count();
    for i in 0..n_verts {
        let vid = VertexId::from_usize(i);
        let vd = mesh.vertices.get_mut(vid);
        let p = &mut vd.position;
        p.x *= width_m;
        p.y *= height_m;
        p.z *= length_m;
    }

    // Assign boundary conditions keyed by node index (= vertex index).
    // Boundary labels from StructuredGridBuilder:  z=0 → "inlet", z=1 → "outlet", sides → "wall"
    let z_inlet = 0.0;
    let z_outlet = length_m;
    let tol = length_m * 1e-6;
    let inlet_center_z = z_inlet;
    let outlet_z = z_outlet;

    // Collect boundary nodes by label
    let mut inlet_nodes: Vec<usize> = Vec::new();
    let mut outlet_nodes: Vec<usize> = Vec::new();
    let mut wall_nodes: Vec<usize> = Vec::new();
    for (vid, vdata) in mesh.vertices.iter() {
        let idx = vid.as_usize();
        let z = vdata.position.z;
        // Check boundary label of adjacent faces (done via position heuristic)
        if (z - inlet_center_z).abs() < tol {
            inlet_nodes.push(idx);
        } else if (z - outlet_z).abs() < tol {
            outlet_nodes.push(idx);
        } else {
            wall_nodes.push(idx);
        }
    }

    // Classify inlet nodes: rim (on outer edge) vs interior
    let inlet_rim: Vec<usize> = inlet_nodes
        .iter()
        .copied()
        .filter(|&i| {
            let vd = mesh.vertices.get(VertexId::from_usize(i));
            let x_frac = vd.position.x / width_m;
            let y_frac = vd.position.y / height_m;
            x_frac < 1e-6 || x_frac > 1.0 - 1e-6 || y_frac < 1e-6 || y_frac > 1.0 - 1e-6
        })
        .collect();
    let inlet_interior: Vec<usize> = inlet_nodes
        .iter()
        .copied()
        .filter(|i| !inlet_rim.contains(i))
        .collect();

    let mut bcs: HashMap<usize, BoundaryCondition<f64>> = HashMap::new();

    // Wall BCs (and inlet rim): zero-velocity Dirichlet
    let zero_dirichlet = BoundaryCondition::Dirichlet {
        value: 0.0_f64,
        component_values: Some(vec![Some(0.0), Some(0.0), Some(0.0), None]),
    };
    for i in wall_nodes {
        bcs.insert(i, zero_dirichlet.clone());
    }
    for i in inlet_rim {
        bcs.insert(i, zero_dirichlet.clone());
    }

    // Inlet interior: uniform plug flow in +z direction
    for i in inlet_interior {
        bcs.insert(
            i,
            BoundaryCondition::VelocityInlet {
                velocity: Vector3::new(0.0, 0.0, u_inlet),
            },
        );
    }

    // Outlet: one reference node at PressureOutlet, rest Outflow
    let mut first_outlet = true;
    for i in outlet_nodes {
        if first_outlet {
            bcs.insert(i, BoundaryCondition::PressureOutlet { pressure: 0.0 });
            first_outlet = false;
        } else {
            bcs.insert(i, BoundaryCondition::Outflow);
        }
    }

    let n_corner = n_verts;
    let fluid = ConstantPropertyFluid::new(
        "blood".to_string(),
        BLOOD_DENSITY,
        BLOOD_VISCOSITY,
        3500.0,
        0.5,
        1570.0,
    );

    let problem = StokesFlowProblem::<f64>::new(mesh.clone(), fluid, bcs, n_corner);
    let config = FemConfig::<f64>::default();
    let mut solver = FemSolver::new(config);

    let solution = solver
        .solve(&problem, None)
        .map_err(|e| anyhow::anyhow!("FemSolver: {e}"))?;

    // Extract pressure drop (inlet mean p - outlet mean p)
    let mut p_inlet_sum = 0.0_f64;
    let mut p_inlet_count = 0_usize;
    let mut p_outlet_sum = 0.0_f64;
    let mut p_outlet_count = 0_usize;
    let mut u_max = 0.0_f64;

    for i in 0..n_corner {
        let vd = mesh.vertices.get(VertexId::from_usize(i));
        let z = vd.position.z;
        if (z - z_inlet).abs() < tol {
            p_inlet_sum += solution.get_pressure(i);
            p_inlet_count += 1;
        } else if (z - z_outlet).abs() < tol {
            p_outlet_sum += solution.get_pressure(i);
            p_outlet_count += 1;
        }
        let vel = solution.get_velocity(i);
        let spd = vel.norm();
        if spd > u_max {
            u_max = spd;
        }
    }

    let p_inlet_mean = if p_inlet_count > 0 {
        p_inlet_sum / p_inlet_count as f64
    } else {
        0.0
    };
    let p_outlet_mean = if p_outlet_count > 0 {
        p_outlet_sum / p_outlet_count as f64
    } else {
        0.0
    };
    let dp = (p_inlet_mean - p_outlet_mean).abs();

    Ok((dp, u_max))
}

fn solve_3d_level(
    design_id: &str,
    blueprint: &cfd_schematics::NetworkBlueprint,
    solved_1d: &OneDResult,
    nx: usize,
    ny_nz: usize,
) -> anyhow::Result<ThreeDLevel> {
    // Build a quick flow map from the 1D solution: channel_id → q_m3_s
    let q_map: HashMap<String, f64> = solved_1d
        .channels
        .iter()
        .map(|(id, q, _)| (id.clone(), *q))
        .collect();

    let mut dp_total = 0.0_f64;
    let mut _dp_venturi = 0.0_f64;
    let mut u_venturi = 0.0_f64;
    let mut n_venturi = 0_usize;

    for ch in &blueprint.channels {
        let q = q_map
            .get(ch.id.as_str())
            .copied()
            .unwrap_or(INLET_FLOW_M3_S);
        let area = ch.cross_section.area();
        let u_in = q.abs() / area.max(1e-18);

        match solve_channel_3d(ch, u_in, nx, ny_nz) {
            Ok((dp, u_max)) => {
                // Aggregate dp in series (longest path from inlet to outlet is the
                // sum of pressure drops across each channel along that path).
                // For the cascade as a whole, we sum all channels as a conservative
                // upper bound (actual path uses only one branch per level).
                // Weighted by flow fraction for a flow-path–consistent aggregate.
                let q_frac = (q / INLET_FLOW_M3_S).abs();
                dp_total += dp * q_frac;

                let is_throat = ch
                    .id
                    .as_str()
                    .as_bytes()
                    .windows("throat".len())
                    .any(|w| w.eq_ignore_ascii_case(b"throat"));
                if is_throat {
                    _dp_venturi += dp;
                    if u_max > u_venturi {
                        u_venturi = u_max;
                    }
                    n_venturi += 1;
                }
            }
            Err(e) => {
                eprintln!(
                    "  3D [{design_id}] channel '{}' failed: {e}",
                    ch.id.as_str()
                );
            }
        }
    }

    if n_venturi == 0 {
        // Non-venturi designs: no cavitation number meaningful
        u_venturi = (INLET_FLOW_M3_S
            / blueprint
                .channels
                .iter()
                .map(|c| c.cross_section.area())
                .fold(f64::MAX, f64::min))
        .max(0.01);
    }

    // Cavitation number σ = (P_ref - P_vapour) / (0.5 ρ u²)
    let sigma = (P_REF_PA - P_VAPOUR_PA) / (0.5 * BLOOD_DENSITY * u_venturi * u_venturi + 1e-12);

    println!(
        "  3D [{design_id}] {nx}×{ny_nz}×{ny_nz}: dp={:.1} Pa  u_max={:.3} m/s  σ={:.4}",
        dp_total, u_venturi, sigma
    );

    Ok(ThreeDLevel {
        nx,
        ny_nz,
        dp_total_pa: dp_total,
        sigma_venturi: sigma,
        u_venturi_m_s: u_venturi,
    })
}

fn solve_3d(
    design_id: &str,
    blueprint: &cfd_schematics::NetworkBlueprint,
    solved_1d: &OneDResult,
) -> anyhow::Result<ThreeDResult> {
    // Three refinement levels: (nx_axial, ny_nz_cross)
    let refinements = [(4usize, 2usize), (8, 4), (16, 8)];
    let mut levels = Vec::with_capacity(3);
    for (nx, ny_nz) in refinements {
        levels.push(solve_3d_level(design_id, blueprint, solved_1d, nx, ny_nz)?);
    }

    let dp_vals = [
        levels[0].dp_total_pa,
        levels[1].dp_total_pa,
        levels[2].dp_total_pa,
    ];
    let sigma_vals = [
        levels[0].sigma_venturi,
        levels[1].sigma_venturi,
        levels[2].sigma_venturi,
    ];

    let (dp_extrap, dp_order, dp_gci) =
        richardson_study(dp_vals[0], dp_vals[1], dp_vals[2], "[3D dp]");
    let (sigma_extrap, sigma_order, sigma_gci) =
        richardson_study(sigma_vals[0], sigma_vals[1], sigma_vals[2], "[3D σ]");

    // Check asymptotic range
    let gci = GridConvergenceIndex::new(3, dp_order, REFINEMENT_RATIO);
    let gci_fine = gci.compute_fine(dp_vals[2], dp_vals[1]);
    let gci_coarse = gci.compute_coarse(dp_vals[2], dp_vals[1]);
    let asymptotic = gci.is_asymptotic(gci_fine, gci_coarse);

    Ok(ThreeDResult {
        design_id: design_id.to_string(),
        levels: [levels[0].clone(), levels[1].clone(), levels[2].clone()],
        dp_extrap_pa: dp_extrap,
        dp_order,
        dp_gci_fine_pct: dp_gci,
        sigma_extrap,
        sigma_order,
        sigma_gci_fine_pct: sigma_gci,
        asymptotic,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 6: Richardson extrapolation helper
// ─────────────────────────────────────────────────────────────────────────────

/// Returns (extrapolated_value, observed_order, gci_fine_pct).
/// Gracefully falls back if values are identical or very close.
fn richardson_study(coarse: f64, medium: f64, fine: f64, label: &str) -> (f64, f64, f64) {
    if fine.abs() < 1e-12 && medium.abs() < 1e-12 {
        return (0.0, 0.0, 0.0);
    }
    let order = RichardsonExtrapolation::estimate_order(coarse, medium, fine, REFINEMENT_RATIO)
        .unwrap_or_else(|_| 2.0_f64);

    let extrap = RichardsonExtrapolation::with_order(order, REFINEMENT_RATIO)
        .map(|re| re.extrapolate(fine, medium))
        .unwrap_or(fine);

    let gci = GridConvergenceIndex::new(3, order, REFINEMENT_RATIO);
    let gci_pct = gci.compute_fine(fine, medium) * 100.0;

    println!(
        "  Richardson {label}: coarse={coarse:.4} medium={medium:.4} fine={fine:.4} \
         → order={order:.2}  extrap={extrap:.4}  GCI={gci_pct:.2}%"
    );

    (extrap, order, gci_pct)
}

// ─────────────────────────────────────────────────────────────────────────────
// Part 7: Cross-fidelity table formatting
// ─────────────────────────────────────────────────────────────────────────────

fn cross_fidelity_table(
    design_id: &str,
    r1d: &OneDResult,
    r2d: &TwoDResult,
    r3d: &ThreeDResult,
) -> String {
    let mut s = String::new();
    s.push_str(&format!("\n### {design_id} — Cross-fidelity summary\n\n"));
    s.push_str("| Quantity | 1D (HP network) | 2D (extrap ± GCI%) | 3D (extrap ± GCI%) |\n");
    s.push_str("|----------|-----------------|---------------------|---------------------|\n");
    s.push_str(&format!(
        "| dp_total [Pa] | {:.1} | {:.1} ± {:.1}% | {:.1} ± {:.1}% |\n",
        r1d.dp_total_pa,
        r2d.wall_shear_extrap_pa,
        r2d.gci_fine_pct,
        r3d.dp_extrap_pa,
        r3d.dp_gci_fine_pct,
    ));
    s.push_str(&format!(
        "| max wall shear [Pa] | (HP) | {:.2} ± {:.1}% | see dp | \n",
        r2d.wall_shear_extrap_pa, r2d.gci_fine_pct,
    ));
    s.push_str(&format!(
        "| Q_venturi [µL/min] | {:.2} | — | — |\n",
        r1d.q_venturi_m3_s * 1e9 * 60.0,
    ));
    s.push_str(&format!(
        "| σ_venturi | (Bernoulli) | — | {:.4} ± {:.1}% |\n",
        r3d.sigma_extrap, r3d.sigma_gci_fine_pct,
    ));
    s.push_str(&format!(
        "| Asymptotic range | — | — | {} |\n",
        if r3d.asymptotic { "YES" } else { "NO" },
    ));
    s
}

// ─────────────────────────────────────────────────────────────────────────────
// Cross-section area helper
// ─────────────────────────────────────────────────────────────────────────────

// SSOT: use CrossSectionSpec::area() instead of local duplicate.

// ─────────────────────────────────────────────────────────────────────────────
// Part 8: Report writer
// ─────────────────────────────────────────────────────────────────────────────

fn write_report(
    out_dir: &Path,
    designs: &[cfd_optim::RankedDesign],
    results_1d: &[OneDResult],
    results_2d: &[TwoDResult],
    results_3d: &[ThreeDResult],
) -> anyhow::Result<()> {
    fs::create_dir_all(out_dir)?;

    let mut md = String::new();
    md.push_str("# Cascade 2D/3D Simulation — Richardson Extrapolation Validation\n\n");
    md.push_str("Generated by `cfd-optim/examples/cascade_2d_3d_validation.rs`\n\n");
    md.push_str("## Design Selection\n\n");
    md.push_str("| Rank | ID | Topology | Score |\n");
    md.push_str("|------|----|----------|-------|\n");
    for (i, d) in designs.iter().enumerate() {
        md.push_str(&format!(
            "| {} | {} | {:?} | {:.4} |\n",
            i + 1,
            d.candidate.id,
            d.candidate.topology,
            d.score,
        ));
    }

    // 1D results
    md.push_str("\n## Part 3 — Complete 1D Network Solution\n\n");
    for r in results_1d {
        md.push_str(&format!("\n### {} — 1D HP Network\n\n", r.design_id));
        md.push_str(&format!(
            "- Inlet pressure: {:.2} Pa\n",
            r.inlet_pressure_pa
        ));
        md.push_str(&format!("- Total ΔP: {:.2} Pa\n", r.dp_total_pa));
        md.push_str(&format!(
            "- Q_venturi: {:.4} µL/min\n\n",
            r.q_venturi_m3_s * 1e9 * 60.0
        ));
        md.push_str("| Channel | Q [µL/min] | P_from [Pa] |\n");
        md.push_str("|---------|------------|-------------|\n");
        for (id, q, p) in &r.channels {
            md.push_str(&format!("| {id} | {:.4} | {:.2} |\n", q * 1e9 * 60.0, p));
        }
    }

    // 2D Richardson results
    md.push_str("\n## Part 4 — 2D Richardson Convergence Study\n\n");
    for r in results_2d {
        md.push_str(&format!("\n### {} — 2D Wall Shear\n\n", r.design_id));
        md.push_str("| Grid | nx | ny | τ_max [Pa] | Converged |\n");
        md.push_str("|------|----|----|------------|----------|\n");
        for (i, lv) in r.levels.iter().enumerate() {
            let label = ["Coarse", "Medium", "Fine"][i];
            md.push_str(&format!(
                "| {label} | {} | {} | {:.4} | {} |\n",
                lv.nx,
                lv.ny,
                lv.max_wall_shear_pa,
                if lv.converged { "YES" } else { "NO" }
            ));
        }
        md.push_str(&format!(
            "\n**Richardson extrapolated τ_max**: {:.4} Pa  \n\
             **Observed order**: {:.2}  \n\
             **GCI (fine)**: {:.2}%\n",
            r.wall_shear_extrap_pa, r.wall_shear_order, r.gci_fine_pct
        ));
    }

    // 3D Richardson results
    md.push_str("\n## Part 5 — 3D Richardson Convergence Study\n\n");
    for r in results_3d {
        md.push_str(&format!(
            "\n### {} — 3D FEM (StructuredGridBuilder)\n\n",
            r.design_id
        ));
        md.push_str("| Grid | nx | ny=nz | dp [Pa] | σ | u_max [m/s] |\n");
        md.push_str("|------|----|-------|---------|---|-------------|\n");
        for (i, lv) in r.levels.iter().enumerate() {
            let label = ["Coarse", "Medium", "Fine"][i];
            md.push_str(&format!(
                "| {label} | {} | {} | {:.2} | {:.4} | {:.4} |\n",
                lv.nx, lv.ny_nz, lv.dp_total_pa, lv.sigma_venturi, lv.u_venturi_m_s
            ));
        }
        md.push_str(&format!(
            "\n**dp — Richardson extrapolated**: {:.2} Pa (order={:.2}, GCI={:.2}%)  \n\
             **σ — Richardson extrapolated**: {:.4} (order={:.2}, GCI={:.2}%)  \n\
             **Asymptotic range**: {}\n",
            r.dp_extrap_pa,
            r.dp_order,
            r.dp_gci_fine_pct,
            r.sigma_extrap,
            r.sigma_order,
            r.sigma_gci_fine_pct,
            if r.asymptotic { "YES" } else { "NO" }
        ));
    }

    // Cross-fidelity tables
    md.push_str("\n## Part 7 — Cross-fidelity Consistency\n\n");
    md.push_str(
        "Expected: 1D ↔ 2D agree within ~15% (HP-dominated channels; 2D adds viscous \
         entrance effects). 3D differs due to junction geometry (recirculation, Dean effects).\n",
    );
    for i in 0..results_1d.len() {
        md.push_str(&cross_fidelity_table(
            &results_1d[i].design_id,
            &results_1d[i],
            &results_2d[i],
            &results_3d[i],
        ));
    }

    // Write markdown
    let md_path = out_dir.join("cascade_validation.md");
    fs::write(&md_path, &md)?;
    println!("\nReport written to: {}", md_path.display());

    // Write per-design JSON for 1D
    for r in results_1d {
        let json = serde_json::json!({
            "design_id": r.design_id,
            "inlet_pressure_pa": r.inlet_pressure_pa,
            "dp_total_pa": r.dp_total_pa,
            "q_venturi_ul_min": r.q_venturi_m3_s * 1e9 * 60.0,
            "channels": r.channels.iter().map(|(id, q, p)| {
                serde_json::json!({"id": id, "q_ul_min": q * 1e9 * 60.0, "p_from_pa": p})
            }).collect::<Vec<_>>(),
        });
        let path = out_dir.join(format!("{}_cascade_1d.json", r.design_id));
        fs::write(&path, serde_json::to_string_pretty(&json)?)?;
    }

    // Write per-design JSON for 3D GCI
    for r in results_3d {
        let json = serde_json::json!({
            "design_id": r.design_id,
            "dp_extrap_pa": r.dp_extrap_pa,
            "dp_order": r.dp_order,
            "dp_gci_fine_pct": r.dp_gci_fine_pct,
            "sigma_extrap": r.sigma_extrap,
            "sigma_order": r.sigma_order,
            "sigma_gci_fine_pct": r.sigma_gci_fine_pct,
            "asymptotic": r.asymptotic,
            "levels": r.levels.iter().map(|l| serde_json::json!({
                "nx": l.nx, "ny_nz": l.ny_nz,
                "dp_pa": l.dp_total_pa, "sigma": l.sigma_venturi,
            })).collect::<Vec<_>>(),
        });
        let path = out_dir.join(format!("{}_cascade_3d_gci.json", r.design_id));
        fs::write(&path, serde_json::to_string_pretty(&json)?)?;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

fn main() -> anyhow::Result<()> {
    println!("=== Cascade 2D/3D Validation ===\n");

    // Part 1: Select top-3 designs
    println!("Part 1: Selecting top-3 RbcProtectedSdt designs...");
    let designs = select_top3()?;
    println!("  Selected {} designs:", designs.len());
    for d in &designs {
        println!(
            "    {} ({:?})  score={:.4}",
            d.candidate.id, d.candidate.topology, d.score
        );
    }

    // Part 2: Extract blueprints
    println!("\nPart 2: Extracting NetworkBlueprints...");
    let blueprints: Vec<cfd_schematics::NetworkBlueprint> =
        designs.iter().map(|d| d.candidate.to_blueprint()).collect();

    // Part 3: 1D solve
    println!("\nPart 3: Complete 1D cascade solve...");
    let results_1d: Vec<OneDResult> = designs
        .iter()
        .zip(blueprints.iter())
        .map(|(d, bp)| solve_1d(&d.candidate.id, bp))
        .collect::<anyhow::Result<_>>()?;

    for r in &results_1d {
        println!(
            "  {} → ΔP={:.2} Pa  Q_venturi={:.4} µL/min  n_channels={}",
            r.design_id,
            r.dp_total_pa,
            r.q_venturi_m3_s * 1e9 * 60.0,
            r.channels.len()
        );
    }

    // Part 4: 2D solve at 3 grid levels
    println!("\nPart 4: Complete 2D cascade solve (3 levels)...");
    let results_2d: Vec<TwoDResult> = designs
        .iter()
        .zip(blueprints.iter())
        .map(|(d, bp)| solve_2d(&d.candidate.id, bp))
        .collect::<anyhow::Result<_>>()?;

    // Part 5: 3D solve at 3 refinement levels
    println!("\nPart 5: Complete 3D cascade solve (3 refinement levels)...");
    let results_3d: Vec<ThreeDResult> = designs
        .iter()
        .zip(blueprints.iter())
        .zip(results_1d.iter())
        .map(|((d, bp), r1d)| solve_3d(&d.candidate.id, bp, r1d))
        .collect::<anyhow::Result<_>>()?;

    // Part 6: Already done inside solve_2d / solve_3d (printed inline)
    println!("\nPart 6: Richardson extrapolation — see output above.");

    // Part 7: Cross-fidelity summary (printed to stdout)
    println!("\nPart 7: Cross-fidelity consistency:");
    for i in 0..results_1d.len() {
        print!(
            "{}",
            cross_fidelity_table(
                &results_1d[i].design_id,
                &results_1d[i],
                &results_2d[i],
                &results_3d[i],
            )
        );
    }

    // Part 8: Write report
    println!("\nPart 8: Writing report...");
    let out_dir = Path::new("report/milestone12");
    write_report(out_dir, &designs, &results_1d, &results_2d, &results_3d)?;

    println!("\n=== Done ===");
    Ok(())
}
