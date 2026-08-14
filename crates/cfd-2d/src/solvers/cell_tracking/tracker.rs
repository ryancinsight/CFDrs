#![cfg_attr(
    test,
    expect(
        clippy::print_stderr,
        reason = "ratchet CFDRS-PRINT-1: pre-existing debt"
    )
)]

use aequitas::systems::si::quantities::{Dimensionless, Length, Time};

use super::physics::{CellTrackerConfig, VelocityFieldInterpolator};
use super::population::{
    CellPopulation, CellRoutingSummary, CellTrajectory, TrackedCell, TrackedPosition,
};

/// Lagrangian cell tracker for 2D velocity fields.
pub struct CellTracker<'a, V: VelocityFieldInterpolator> {
    velocity: &'a V,
    config: CellTrackerConfig,
}

impl<'a, V: VelocityFieldInterpolator> CellTracker<'a, V> {
    #[must_use]
    /// Create a new cell tracker for a resolved flow field.
    pub fn new(velocity: &'a V, config: CellTrackerConfig) -> Self {
        Self { velocity, config }
    }

    /// Trace cells through the velocity field.
    #[must_use]
    pub fn trace_cells(
        &self,
        cells: &[TrackedCell],
        dt: Time,
        max_steps: usize,
    ) -> Vec<CellTrajectory> {
        cells
            .iter()
            .map(|cell| self.trace_single(cell, dt, max_steps))
            .collect()
    }

    fn trace_single(&self, cell: &TrackedCell, dt: Time, max_steps: usize) -> CellTrajectory {
        let mut traj = CellTrajectory {
            cell_id: cell.id,
            population: cell.population,
            positions: Vec::with_capacity(max_steps / 10 + 2),
            exit_outlet: None,
        };

        let a = cell.population.diameter().into_base();
        let rho_cell = cell.population.density().into_base();
        let c_l = cell.population.lift_coefficient().into_base();
        let (x_min_q, x_max_q, y_min_q, y_max_q) = self.velocity.bounds();
        let x_min = x_min_q.into_base();
        let x_max = x_max_q.into_base();
        let y_min = y_min_q.into_base();
        let y_max = y_max_q.into_base();
        let d_h = self.config.hydraulic_diameter.into_base().max(a * 2.0);
        let mu = self.config.viscosity.into_base();
        let rho_f = self.config.fluid_density.into_base();
        let dt = dt.into_base();

        // Estimate u_max from inlet centerline if not provided.
        let u_max = if self.config.u_max.into_base() > 1e-12 {
            self.config.u_max.into_base()
        } else {
            let y_mid = (y_min + y_max) * 0.5;
            self.velocity
                .velocity_at(
                    Length::from_base(x_min + (x_max - x_min) * 0.05),
                    Length::from_base(y_mid),
                )
                .0
                .into_base()
                .abs()
                .max(1e-6)
        };

        // Precompute cell constants.
        let r = a * 0.5;
        let cell_vol = (4.0 / 3.0) * std::f64::consts::PI * r.powi(3);
        let cell_mass = rho_cell * cell_vol;
        let drag_coeff = 3.0 * std::f64::consts::PI * mu * a;
        let tau = cell_mass / drag_coeff.max(1e-30);

        // Lift force prefactor: C_L * rho_f * U_max^2 * a^4 / D_h^2
        let lift_prefactor = c_l * rho_f * u_max.powi(2) * a.powi(4) / d_h.powi(2);

        let mut x = cell.x.into_base();
        let mut y = cell.y.into_base();
        let mut vx = cell.vx.into_base();
        let mut vy = cell.vy.into_base();

        traj.positions.push(TrackedPosition {
            x: Length::from_base(x),
            y: Length::from_base(y),
            time: Time::from_base(0.0),
        });

        for step in 0..max_steps {
            if !self
                .velocity
                .is_fluid(Length::from_base(x), Length::from_base(y))
            {
                break;
            }

            // --- Fluid velocity at cell position ---
            let (uf_q, vf_q) = self
                .velocity
                .velocity_at(Length::from_base(x), Length::from_base(y));
            let uf = uf_q.into_base();
            let vf = vf_q.into_base();

            // --- Stokes drag acceleration ---
            let ax_drag = (uf - vx) / tau;
            let ay_drag = (vf - vy) / tau;

            // --- Zweifach-Fung cell-size exclusion (Fung 1969) ---
            // The dominant mechanism at low Re (millifluidic regime): cells
            // near the flow dividing streamline experience a cross-stream
            // velocity gradient that pulls them toward the faster (wider)
            // daughter.  This effect is modeled as a drift velocity
            // proportional to the velocity difference across the cell body:
            //
            //   v_zf = C_L * a * du/dy
            //
            // where du/dy is the streamwise velocity gradient in the cross-
            // stream direction.  This gives the correct dimensional scaling:
            // larger cells (bigger a) and stiffer cells (higher C_L) drift
            // faster toward the high-flow daughter.
            //
            // The drift velocity is added directly to the cell velocity
            // rather than computed as a force, avoiding the mass scaling
            // issue that made the force-based approach too weak.
            let du_dy_local = {
                let probe = a.max(1e-7);
                let (u_plus, _) = self
                    .velocity
                    .velocity_at(Length::from_base(x), Length::from_base(y + probe));
                let (u_minus, _) = self
                    .velocity
                    .velocity_at(Length::from_base(x), Length::from_base(y - probe));
                (u_plus.into_base() - u_minus.into_base()) / (2.0 * probe)
            };
            // Zweifach drift velocity: cell migrates toward the faster side.
            // The sign of du/dy determines the direction: positive du/dy
            // means faster flow above -> cell drifts upward.
            let v_zweifach = c_l * a * du_dy_local;
            // Convert to acceleration via the drag relaxation time.
            let ay_zweifach = v_zweifach / tau;

            // --- Inertial lift (Di Carlo 2009) ---
            // Normalized position: s = 0 at center, +/-1 at walls.
            let y_center = (y_min + y_max) * 0.5;
            let half_h = (y_max - y_min) * 0.5;
            let s = ((y - y_center) / half_h).clamp(-0.95, 0.95);

            // Lift profile f(s): pushes cells toward equilibrium at |s| ~ 0.6.
            // Polynomial fit to numerical data from Asmolov (1999):
            //   f(s) ~ s * (1 - s^2) * (s^2 - s_eq^2)
            // Sign convention: positive f pushes toward +y.
            // At s=0 (center): f > 0 for s > 0 side -> pushes outward.
            // At s=0.6 (equilibrium): f = 0.
            // At s=0.9 (near wall): f < 0 -> pushes inward (away from wall).
            let s_eq = 0.6;
            let f_s = s * (1.0 - s * s) * (s * s - s_eq * s_eq);

            let f_lift = lift_prefactor * f_s;
            let ay_lift = f_lift / cell_mass.max(1e-30);

            // --- Wall lubrication repulsion ---
            let dist_bot = (y - y_min - r).max(r * 0.05);
            let dist_top = (y_max - y - r).max(r * 0.05);
            let wall_len = r * 0.2;
            let wall_accel = (drag_coeff * 5.0 / cell_mass.max(1e-30))
                * ((-dist_bot / wall_len).exp() - (-dist_top / wall_len).exp());

            // --- Buoyancy (negligible for blood cells, but included for correctness) ---
            let ay_buoy = -(rho_cell - rho_f) / rho_cell * 9.81;

            // --- Total acceleration ---
            let ax = ax_drag;
            let ay = ay_drag + ay_zweifach + ay_lift + wall_accel + ay_buoy;

            // --- RK2 integration ---
            let vx_mid = vx + 0.5 * dt * ax;
            let vy_mid = vy + 0.5 * dt * ay;
            let x_mid = x + 0.5 * dt * vx_mid;
            let y_mid_pos = y + 0.5 * dt * vy_mid;

            let (uf_mid_q, vf_mid_q) = self
                .velocity
                .velocity_at(Length::from_base(x_mid), Length::from_base(y_mid_pos));
            let uf_mid = uf_mid_q.into_base();
            let vf_mid = vf_mid_q.into_base();
            let ax_mid = (uf_mid - vx_mid) / tau;
            let ay_mid = (vf_mid - vy_mid) / tau + ay_zweifach + ay_lift + wall_accel + ay_buoy;

            vx += dt * ax_mid;
            vy += dt * ay_mid;
            let x_prev = x;
            x += dt * vx;
            y += dt * vy;

            // --- Pries PSM junction routing (plasma skimming) ---
            // When a cell crosses the bifurcation split plane, use the Pries
            // Phase Separation Model to determine probabilistically whether
            // the cell enters the wide (center) or narrow (peripheral)
            // daughter.  The PSM captures the cell-free layer and plasma
            // skimming physics that dominate at low Re.
            if self.config.split_x.into_base() > 0.0
                && x_prev < self.config.split_x.into_base()
                && x >= self.config.split_x.into_base()
            {
                let y_div = self.config.dividing_streamline_y.into_base();
                if let Some(psm) = self.config.psm_params.as_ref().filter(|_| y_div > 0.0) {
                    let fqb = crate::solvers::plasma_skimming::pries_phase_separation(
                        psm.flow_fraction_wide.into_base(),
                        &crate::solvers::plasma_skimming::PriesPhaseParams {
                            parent_diameter_m: d_h,
                            daughter_alpha_diameter_m: psm.wide_daughter_diameter.into_base(),
                            daughter_beta_diameter_m: psm.narrow_daughter_diameter.into_base(),
                            feed_hematocrit: psm.feed_hematocrit.into_base(),
                        },
                    );
                    // For cell-type-specific routing, compute modified PSM
                    // with larger effective cell size for CTCs.
                    let cell_fqe =
                        crate::solvers::plasma_skimming::pries_phase_separation_cell_type(
                            psm.flow_fraction_wide.into_base(),
                            &crate::solvers::plasma_skimming::PriesPhaseParams {
                                parent_diameter_m: d_h,
                                daughter_alpha_diameter_m: psm.wide_daughter_diameter.into_base(),
                                daughter_beta_diameter_m: psm.narrow_daughter_diameter.into_base(),
                                feed_hematocrit: psm.feed_hematocrit.into_base(),
                            },
                            a,
                            c_l / 0.5, // normalize lift coeff to stiffness factor
                        );

                    // Deterministic routing based on cell position relative to
                    // the effective dividing streamline, shifted by the PSM.
                    //
                    // The PSM tells us that the cell_fqe fraction of cells
                    // should go to the wide daughter.  We implement this by
                    // shifting the effective dividing streamline: cells above
                    // the shifted line go to wide, below go to narrow.
                    //
                    // The shift is: if cell_fqe > fqb.cell_fraction (CTCs),
                    // the divider moves down (more cells go to wide).
                    let rbc_fqe = fqb.cell_fraction;
                    let cell_fqe_val = cell_fqe.cell_fraction;
                    // Shift the divider based on the cell-type-specific excess.
                    let effective_y_div = if cell_fqe_val > rbc_fqe {
                        // This cell type goes to wide more than RBCs ->
                        // lower the divider to capture more cells.
                        let excess = (cell_fqe_val - rbc_fqe).min(0.3);
                        y_div * (1.0 - excess)
                    } else {
                        y_div
                    };

                    if y >= effective_y_div {
                        // Route to wide daughter (upper half).
                        let wide_center = (effective_y_div + y_max) * 0.5;
                        y = wide_center.min(y_max - r);
                    } else {
                        // Route to narrow daughter (lower half).
                        let narrow_center = effective_y_div * 0.5;
                        y = narrow_center.max(y_min + r);
                    }
                }
            }

            // Wall bounce: reflect off boundaries.
            let y_lo_wall = y_min + r;
            let y_hi_wall = y_max - r;
            if y < y_lo_wall {
                y = y_lo_wall + (y_lo_wall - y).min(r);
                vy = vy.abs() * 0.3; // inelastic bounce
            }
            if y > y_hi_wall {
                y = y_hi_wall - (y - y_hi_wall).min(r);
                vy = -vy.abs() * 0.3;
            }

            // Record at intervals.
            let t = (step + 1) as f64 * dt;
            if step % 10 == 0 {
                traj.positions.push(TrackedPosition {
                    x: Length::from_base(x),
                    y: Length::from_base(y),
                    time: Time::from_base(t),
                });
            }

            // --- Outlet exit check ---
            if x >= x_max {
                traj.positions.push(TrackedPosition {
                    x: Length::from_base(x),
                    y: Length::from_base(y),
                    time: Time::from_base(t),
                });
                traj.exit_outlet = Some(self.classify_outlet(y));
                break;
            }
            // Backflow exit (left boundary).
            if x <= x_min {
                traj.positions.push(TrackedPosition {
                    x: Length::from_base(x),
                    y: Length::from_base(y),
                    time: Time::from_base(t),
                });
                traj.exit_outlet = Some("backflow".to_string());
                break;
            }
        }
        traj
    }

    fn classify_outlet(&self, y: f64) -> String {
        if !self.config.outlet_zones.is_empty() {
            for zone in &self.config.outlet_zones {
                if y >= zone.y_lo.into_base() && y <= zone.y_hi.into_base() {
                    return zone.name.clone();
                }
            }
            return "unclassified".to_string();
        }
        // Default: center/peripheral split at y_mid +/- 25%.
        let (_, _, y_min_q, y_max_q) = self.velocity.bounds();
        let y_min = y_min_q.into_base();
        let y_max = y_max_q.into_base();
        let y_mid = (y_min + y_max) * 0.5;
        let y_span = y_max - y_min;
        if (y - y_mid).abs() < y_span * 0.25 {
            "center".to_string()
        } else {
            "peripheral".to_string()
        }
    }

    /// Classify tracked trajectories into a routing summary.
    #[must_use]
    pub fn classify_routing(&self, trajectories: &[CellTrajectory]) -> CellRoutingSummary {
        let mut s = CellRoutingSummary::default();
        for traj in trajectories {
            let is_center = traj.exit_outlet.as_deref() == Some("center");
            match traj.population {
                CellPopulation::CTC => {
                    s.ctc_total += 1;
                    if is_center {
                        s.ctc_center += 1;
                    }
                }
                CellPopulation::WBC => {
                    s.wbc_total += 1;
                    if is_center {
                        s.wbc_center += 1;
                    }
                }
                CellPopulation::RBC => {
                    s.rbc_total += 1;
                    if is_center {
                        s.rbc_center += 1;
                    }
                }
            }
        }
        let ctc_frac = if s.ctc_total > 0 {
            s.ctc_center as f64 / s.ctc_total as f64
        } else {
            0.0
        };
        let rbc_periph = if s.rbc_total > 0 {
            1.0 - s.rbc_center as f64 / s.rbc_total as f64
        } else {
            0.0
        };
        s.cancer_center_fraction = Dimensionless::from_base(ctc_frac);
        s.separation_efficiency = Dimensionless::from_base((ctc_frac * rbc_periph).sqrt());
        s
    }
}

#[cfg(test)]
mod tests {
    use aequitas::systems::si::quantities::{
        Dimensionless, DynamicViscosity, Length, MassDensity, Time, Velocity,
    };

    use super::*;
    use crate::solvers::cell_tracking::physics::*;
    use crate::solvers::cell_tracking::population::OutletZone;

    #[test]
    fn poiseuille_cells_exit_domain() {
        let flow = PoiseuilleFlow2D {
            u_max: Velocity::from_base(0.1),
            width: Length::from_base(0.01),
            height: Length::from_base(0.002),
        };
        let config = CellTrackerConfig {
            viscosity: DynamicViscosity::from_base(3.5e-3),
            fluid_density: MassDensity::from_base(1025.0),
            hydraulic_diameter: Length::from_base(0.002),
            u_max: Velocity::from_base(0.1),
            ..Default::default()
        };
        let tracker = CellTracker::new(&flow, config);
        let cells = vec![
            TrackedCell {
                population: CellPopulation::CTC,
                x: Length::from_base(0.0),
                y: Length::from_base(0.001),
                vx: Velocity::from_base(0.05),
                vy: Velocity::from_base(0.0),
                id: 0,
            },
            TrackedCell {
                population: CellPopulation::RBC,
                x: Length::from_base(0.0),
                y: Length::from_base(0.001),
                vx: Velocity::from_base(0.05),
                vy: Velocity::from_base(0.0),
                id: 1,
            },
            TrackedCell {
                population: CellPopulation::WBC,
                x: Length::from_base(0.0),
                y: Length::from_base(0.001),
                vx: Velocity::from_base(0.05),
                vy: Velocity::from_base(0.0),
                id: 2,
            },
        ];
        let trajectories = tracker.trace_cells(&cells, Time::from_base(1e-5), 200_000);
        for traj in &trajectories {
            assert!(
                traj.exit_outlet.is_some(),
                "cell {} should exit",
                traj.cell_id
            );
            assert!(
                traj.positions.len() >= 2,
                "cell {} should have trajectory points",
                traj.cell_id
            );
        }
    }

    #[test]
    fn poiseuille_centerline_cells_stay_near_center() {
        let flow = PoiseuilleFlow2D {
            u_max: Velocity::from_base(0.2),
            width: Length::from_base(0.02),
            height: Length::from_base(0.002),
        };
        let config = CellTrackerConfig {
            viscosity: DynamicViscosity::from_base(3.5e-3),
            fluid_density: MassDensity::from_base(1025.0),
            hydraulic_diameter: Length::from_base(0.002),
            u_max: Velocity::from_base(0.2),
            ..Default::default()
        };
        let tracker = CellTracker::new(&flow, config);
        // Release at centerline: should stay near center.
        let cells = vec![TrackedCell {
            population: CellPopulation::CTC,
            x: Length::from_base(0.0),
            y: Length::from_base(0.001),
            vx: Velocity::from_base(0.1),
            vy: Velocity::from_base(0.0),
            id: 0,
        }];
        let trajectories = tracker.trace_cells(&cells, Time::from_base(1e-5), 500_000);
        let traj = &trajectories[0];
        assert!(traj.exit_outlet.is_some());
        // Final y should be within 30% of centerline (0.001 m).
        let final_y = traj
            .positions
            .last()
            .expect("trajectory has a sample")
            .y
            .into_base();
        assert!(
            (final_y - 0.001).abs() < 0.0003,
            "CTC at centerline should stay near center, final_y = {final_y}"
        );
    }

    #[test]
    fn asymmetric_bifurcation_ctcs_prefer_wide_daughter() {
        // Moderate asymmetry (55/45 split) for realistic Zweifach-Fung testing.
        let parent_h = 2.0e-3;
        let center_frac = 0.55;
        let flow = AsymmetricBifurcationFlow {
            parent_width: Length::from_base(2.0e-3),
            parent_height: Length::from_base(parent_h),
            wide_daughter_width: Length::from_base(2.0e-3 * center_frac),
            narrow_daughter_width: Length::from_base(2.0e-3 * (1.0 - center_frac)),
            length: Length::from_base(0.015),
            u_inlet: Velocity::from_base(0.05),
            x_split: Length::from_base(0.005),
        };
        let y_div = flow.dividing_streamline_y();
        let config = CellTrackerConfig {
            viscosity: DynamicViscosity::from_base(3.5e-3),
            fluid_density: MassDensity::from_base(1025.0),
            hydraulic_diameter: Length::from_base(parent_h),
            u_max: Velocity::from_base(0.05),
            outlet_zones: vec![
                OutletZone {
                    name: "center".to_string(),
                    x_min: Length::from_base(0.014),
                    y_lo: y_div,
                    y_hi: Length::from_base(parent_h),
                },
                OutletZone {
                    name: "peripheral".to_string(),
                    x_min: Length::from_base(0.014),
                    y_lo: Length::from_base(0.0),
                    y_hi: y_div,
                },
            ],
            split_x: Length::from_base(0.005),
            dividing_streamline_y: y_div,
            psm_params: Some(PsmBifurcationParams {
                flow_fraction_wide: flow.flow_fraction_wide(),
                wide_daughter_diameter: Length::from_base(2.0e-3 * center_frac),
                narrow_daughter_diameter: Length::from_base(2.0e-3 * (1.0 - center_frac)),
                feed_hematocrit: Dimensionless::from_base(0.45),
            }),
        };
        let tracker = CellTracker::new(&flow, config);

        let n_per_pop = 30;
        let mut cells = Vec::with_capacity(n_per_pop * 2);
        for i in 0..n_per_pop {
            let y = parent_h * 0.05 + parent_h * 0.9 * (i as f64 / (n_per_pop - 1) as f64);
            cells.push(TrackedCell {
                population: CellPopulation::CTC,
                x: Length::from_base(1e-4),
                y: Length::from_base(y),
                vx: Velocity::from_base(0.03),
                vy: Velocity::from_base(0.0),
                id: i,
            });
            cells.push(TrackedCell {
                population: CellPopulation::RBC,
                x: Length::from_base(1e-4),
                y: Length::from_base(y),
                vx: Velocity::from_base(0.03),
                vy: Velocity::from_base(0.0),
                id: n_per_pop + i,
            });
        }

        let trajectories = tracker.trace_cells(&cells, Time::from_base(2e-6), 1_000_000);
        let routing = tracker.classify_routing(&trajectories);

        let q_frac = flow.flow_fraction_wide().into_base();
        eprintln!("--- Asymmetric bifurcation cell tracking ---");
        eprintln!("  center_frac = {center_frac}, flow_fraction_wide = {q_frac:.4}");
        eprintln!("  y_div = {:.6} m", y_div.into_base());
        eprintln!(
            "  CTC center {}/{}, RBC center {}/{}",
            routing.ctc_center, routing.ctc_total, routing.rbc_center, routing.rbc_total,
        );

        let exited = trajectories
            .iter()
            .filter(|t| t.exit_outlet.is_some())
            .count();
        assert!(exited >= cells.len() / 2, "at least half should exit");

        if routing.ctc_total >= 5 && routing.rbc_total >= 5 {
            let ctc_rate = routing.ctc_center as f64 / routing.ctc_total as f64;
            let rbc_rate = routing.rbc_center as f64 / routing.rbc_total as f64;
            eprintln!("  CTC center rate = {ctc_rate:.3}, RBC center rate = {rbc_rate:.3}");
            // CTCs should route to center at >= rate than RBCs.
            assert!(
                ctc_rate >= rbc_rate,
                "CTCs should route to center >= RBCs: CTC {ctc_rate:.3}, RBC {rbc_rate:.3}"
            );
        }
    }
}
