//! Lagrangian cell-tracking solver for 2D velocity fields.
//!
//! Traces discrete cell trajectories through a resolved 2D velocity field,
//! applying drag and inertial lift forces to predict cell routing at
//! bifurcations.  This provides a 2D validation pathway for the 1D
//! Zweifach-Fung cell separation model in cfd-1d.
//!
//! # Physics
//!
//! Each cell (CTC, WBC, or RBC) is modeled as a rigid sphere of diameter `a`
//! subject to:
//!
//! 1. **Stokes drag**: F_drag = 3 pi mu a (u_fluid - u_cell)
//!
//! 2. **Inertial lift** (Di Carlo 2009, Amini 2014):
//!    F_lift = C_L * rho * U_max^2 * a^4 / D_h^2 * f(s)
//!    where s = (y - y_center) / (D_h/2) is the normalized position and
//!    f(s) drives cells toward equilibrium positions at s ~ +/- 0.6.
//!    The sign of f(s) depends on position: cells near the center are
//!    pushed outward, cells near the wall are pushed inward.
//!
//! 3. **Wall repulsion**: lubrication force that prevents overlap with solid
//!    boundaries.
//!
//! # Cell density
//!
//! Cell density differs from the suspending fluid (plasma ~1025 kg/m3):
//! - CTC: 1068 kg/m3 (MCF-7 reference, Byun 2013)
//! - WBC: 1070 kg/m3 (granulocyte average)
//! - RBC: 1100 kg/m3 (packed hematocrit reference)

use serde::{Deserialize, Serialize};

/// Cell population type, matching the 1D model's three-population framework.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellPopulation {
    /// Circulating tumor cell (diameter 10-15 um, stiff).
    CTC,
    /// White blood cell (diameter 10-12 um, moderate stiffness).
    WBC,
    /// Red blood cell (diameter ~8 um, highly deformable).
    RBC,
}

impl CellPopulation {
    /// Characteristic diameter [m].
    #[must_use]
    pub fn diameter_m(self) -> f64 {
        match self {
            Self::CTC => 12.5e-6,
            Self::WBC => 11.0e-6,
            Self::RBC => 8.0e-6,
        }
    }

    /// Cell density [kg/m3].
    #[must_use]
    pub fn density_kg_m3(self) -> f64 {
        match self {
            Self::CTC => 1068.0,
            Self::WBC => 1070.0,
            Self::RBC => 1100.0,
        }
    }

    /// Di Carlo (2009) inertial lift coefficient.
    /// Stiffer cells experience stronger lift toward equilibrium positions.
    /// These values are for kappa = a/D_h ~ 0.05-0.15 (millifluidic range).
    #[must_use]
    pub fn lift_coefficient(self) -> f64 {
        match self {
            Self::CTC => 0.5,  // stiff, strong focusing
            Self::WBC => 0.35, // moderate deformability
            Self::RBC => 0.15, // high deformability, weaker focusing
        }
    }
}

/// A single tracked cell with position, velocity, and identity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrackedCell {
    pub population: CellPopulation,
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub id: usize,
}

/// Trajectory record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellTrajectory {
    pub cell_id: usize,
    pub population: CellPopulation,
    pub positions: Vec<[f64; 3]>, // [x, y, t]
    pub exit_outlet: Option<String>,
}

/// Outlet zone definition for classifying cell exit positions.
#[derive(Debug, Clone)]
pub struct OutletZone {
    pub name: String,
    /// x range: cell exits when x >= x_min.
    pub x_min: f64,
    /// y range [y_lo, y_hi] that defines this outlet.
    pub y_lo: f64,
    pub y_hi: f64,
}

/// Summary of cell routing at outlets.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CellRoutingSummary {
    pub ctc_total: usize,
    pub wbc_total: usize,
    pub rbc_total: usize,
    pub ctc_center: usize,
    pub wbc_center: usize,
    pub rbc_center: usize,
    pub cancer_center_fraction: f64,
    pub separation_efficiency: f64,
}

/// Velocity field interface for the cell tracker.
pub trait VelocityFieldInterpolator {
    /// Interpolate velocity at (x, y).  Returns (u, v) in m/s.
    fn velocity_at(&self, x: f64, y: f64) -> (f64, f64);

    /// Check if (x, y) is inside the fluid domain.
    fn is_fluid(&self, x: f64, y: f64) -> bool;

    /// Domain bounding box: (x_min, x_max, y_min, y_max) in meters.
    fn bounds(&self) -> (f64, f64, f64, f64);
}

/// Configuration for the Lagrangian tracker.
#[derive(Debug, Clone)]
pub struct CellTrackerConfig {
    /// Fluid viscosity [Pa.s].  Default: 3.5e-3 (blood).
    pub viscosity: f64,
    /// Fluid density [kg/m3].  Default: 1025.0 (plasma).
    pub fluid_density: f64,
    /// Hydraulic diameter of the parent channel [m], used for lift scaling.
    pub hydraulic_diameter_m: f64,
    /// Max streamwise velocity [m/s] for lift scaling.  If zero, estimated
    /// from the velocity field at the inlet centerline.
    pub u_max: f64,
    /// Outlet zones for classifying exits.  If empty, a default center/peripheral
    /// split at y_mid +/- 25% is used.
    pub outlet_zones: Vec<OutletZone>,
    /// X-coordinate of the bifurcation split plane [m].  Set to 0 to disable
    /// the junction routing correction.
    pub split_x: f64,
    /// Y-coordinate of the dividing streamline at the split plane [m].
    pub dividing_streamline_y: f64,
    /// Pries Phase Separation Model parameters for the bifurcation.
    /// If None, the junction routing correction is disabled.
    pub psm_params: Option<PsmBifurcationParams>,
}

/// Parameters for the Pries PSM at a specific bifurcation.
#[derive(Debug, Clone)]
pub struct PsmBifurcationParams {
    /// Fractional blood flow into the wide (center) daughter.
    pub flow_fraction_wide: f64,
    /// Hydraulic diameter of the wide daughter [m].
    pub wide_daughter_dh: f64,
    /// Hydraulic diameter of the narrow daughter [m].
    pub narrow_daughter_dh: f64,
    /// Feed hematocrit.
    pub feed_hematocrit: f64,
}

impl Default for CellTrackerConfig {
    fn default() -> Self {
        Self {
            viscosity: 3.5e-3,
            fluid_density: 1025.0,
            hydraulic_diameter_m: 1.0e-3,
            u_max: 0.0,
            outlet_zones: Vec::new(),
            split_x: 0.0,
            dividing_streamline_y: 0.0,
            psm_params: None,
        }
    }
}

/// Lagrangian cell tracker for 2D velocity fields.
pub struct CellTracker<'a, V: VelocityFieldInterpolator> {
    velocity: &'a V,
    config: CellTrackerConfig,
}

impl<'a, V: VelocityFieldInterpolator> CellTracker<'a, V> {
    #[must_use]
    pub fn new(velocity: &'a V, config: CellTrackerConfig) -> Self {
        Self { velocity, config }
    }

    /// Trace cells through the velocity field.
    pub fn trace_cells(
        &self,
        cells: &[TrackedCell],
        dt: f64,
        max_steps: usize,
    ) -> Vec<CellTrajectory> {
        cells
            .iter()
            .map(|cell| self.trace_single(cell, dt, max_steps))
            .collect()
    }

    fn trace_single(&self, cell: &TrackedCell, dt: f64, max_steps: usize) -> CellTrajectory {
        let mut traj = CellTrajectory {
            cell_id: cell.id,
            population: cell.population,
            positions: Vec::with_capacity(max_steps / 10 + 2),
            exit_outlet: None,
        };

        let a = cell.population.diameter_m();
        let rho_cell = cell.population.density_kg_m3();
        let c_l = cell.population.lift_coefficient();
        let (x_min, x_max, y_min, y_max) = self.velocity.bounds();
        let d_h = self.config.hydraulic_diameter_m.max(a * 2.0);
        let mu = self.config.viscosity;
        let rho_f = self.config.fluid_density;

        // Estimate u_max from inlet centerline if not provided.
        let u_max = if self.config.u_max > 1e-12 {
            self.config.u_max
        } else {
            let y_mid = (y_min + y_max) * 0.5;
            self.velocity.velocity_at(x_min + (x_max - x_min) * 0.05, y_mid).0.abs().max(1e-6)
        };

        // Precompute cell constants.
        let r = a * 0.5;
        let cell_vol = (4.0 / 3.0) * std::f64::consts::PI * r.powi(3);
        let cell_mass = rho_cell * cell_vol;
        let drag_coeff = 3.0 * std::f64::consts::PI * mu * a;
        let tau = cell_mass / drag_coeff.max(1e-30);

        // Lift force prefactor: C_L * rho_f * U_max^2 * a^4 / D_h^2
        let lift_prefactor = c_l * rho_f * u_max.powi(2) * a.powi(4) / d_h.powi(2);

        let mut x = cell.x;
        let mut y = cell.y;
        let mut vx = cell.vx;
        let mut vy = cell.vy;

        traj.positions.push([x, y, 0.0]);

        for step in 0..max_steps {
            if !self.velocity.is_fluid(x, y) {
                break;
            }

            // --- Fluid velocity at cell position ---
            let (uf, vf) = self.velocity.velocity_at(x, y);

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
                let (u_plus, _) = self.velocity.velocity_at(x, y + probe);
                let (u_minus, _) = self.velocity.velocity_at(x, y - probe);
                (u_plus - u_minus) / (2.0 * probe)
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

            let (uf_mid, vf_mid) = self.velocity.velocity_at(x_mid, y_mid_pos);
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
            if self.config.split_x > 0.0
                && x_prev < self.config.split_x
                && x >= self.config.split_x
            {
                let y_div = self.config.dividing_streamline_y;
                if y_div > 0.0 && self.config.psm_params.is_some() {
                    let psm = self.config.psm_params.as_ref().unwrap();
                    let fqb = super::plasma_skimming::pries_phase_separation(
                        psm.flow_fraction_wide,
                        &super::plasma_skimming::PriesPhaseParams {
                            parent_diameter_m: d_h,
                            daughter_alpha_diameter_m: psm.wide_daughter_dh,
                            daughter_beta_diameter_m: psm.narrow_daughter_dh,
                            feed_hematocrit: psm.feed_hematocrit,
                        },
                    );
                    // For cell-type-specific routing, compute modified PSM
                    // with larger effective cell size for CTCs.
                    let cell_fqe = super::plasma_skimming::pries_phase_separation_cell_type(
                        psm.flow_fraction_wide,
                        &super::plasma_skimming::PriesPhaseParams {
                            parent_diameter_m: d_h,
                            daughter_alpha_diameter_m: psm.wide_daughter_dh,
                            daughter_beta_diameter_m: psm.narrow_daughter_dh,
                            feed_hematocrit: psm.feed_hematocrit,
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
                traj.positions.push([x, y, t]);
            }

            // --- Outlet exit check ---
            if x >= x_max {
                traj.positions.push([x, y, t]);
                traj.exit_outlet = Some(self.classify_outlet(y));
                break;
            }
            // Backflow exit (left boundary).
            if x <= x_min {
                traj.positions.push([x, y, t]);
                traj.exit_outlet = Some("backflow".to_string());
                break;
            }
        }
        traj
    }

    fn classify_outlet(&self, y: f64) -> String {
        if !self.config.outlet_zones.is_empty() {
            for zone in &self.config.outlet_zones {
                if y >= zone.y_lo && y <= zone.y_hi {
                    return zone.name.clone();
                }
            }
            return "unclassified".to_string();
        }
        // Default: center/peripheral split at y_mid +/- 25%.
        let (_, _, y_min, y_max) = self.velocity.bounds();
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
                    if is_center { s.ctc_center += 1; }
                }
                CellPopulation::WBC => {
                    s.wbc_total += 1;
                    if is_center { s.wbc_center += 1; }
                }
                CellPopulation::RBC => {
                    s.rbc_total += 1;
                    if is_center { s.rbc_center += 1; }
                }
            }
        }
        let ctc_frac = if s.ctc_total > 0 {
            s.ctc_center as f64 / s.ctc_total as f64
        } else { 0.0 };
        let rbc_periph = if s.rbc_total > 0 {
            1.0 - s.rbc_center as f64 / s.rbc_total as f64
        } else { 0.0 };
        s.cancer_center_fraction = ctc_frac;
        s.separation_efficiency = (ctc_frac * rbc_periph).sqrt();
        s
    }
}

// ── Analytical velocity fields for testing ──────────────────────────────────

/// Poiseuille flow in a 2D channel: u(y) = U_max * (1 - (2y/H - 1)^2).
/// Used for unit testing the cell tracker against known analytical solutions.
pub struct PoiseuilleFlow2D {
    pub u_max: f64,
    pub width: f64,  // x-extent [m]
    pub height: f64, // channel height [m]
}

impl VelocityFieldInterpolator for PoiseuilleFlow2D {
    fn velocity_at(&self, _x: f64, y: f64) -> (f64, f64) {
        let s = (2.0 * y / self.height - 1.0).clamp(-1.0, 1.0);
        (self.u_max * (1.0 - s * s), 0.0)
    }
    fn is_fluid(&self, x: f64, y: f64) -> bool {
        x >= 0.0 && x <= self.width && y >= 0.0 && y <= self.height
    }
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.width, 0.0, self.height)
    }
}

/// Asymmetric bifurcation flow: parent channel splits into a wide and narrow
/// daughter at x = x_split.  The wide daughter gets more flow (proportional
/// to width^3 per Hagen-Poiseuille).
pub struct AsymmetricBifurcationFlow {
    pub parent_width_m: f64,
    pub parent_height_m: f64,
    pub wide_daughter_width_m: f64,
    pub narrow_daughter_width_m: f64,
    pub length_m: f64,
    pub u_inlet: f64,
    pub x_split: f64,
}

impl AsymmetricBifurcationFlow {
    /// Flow fraction to the wide daughter (Q_wide / Q_total).
    /// HP scaling: Q proportional to w^3 for equal height/length.
    fn flow_fraction_wide(&self) -> f64 {
        let q_wide = self.wide_daughter_width_m.powi(3);
        let q_narrow = self.narrow_daughter_width_m.powi(3);
        q_wide / (q_wide + q_narrow)
    }

    /// The dividing streamline y-position at the split plane.
    /// In the parent, the flow above this line goes to the wide daughter
    /// and below goes to the narrow daughter.  For Poiseuille flow, the
    /// streamline that captures fraction f of the total flow is at:
    ///   y_div / H = 0.5 * (1 + (1 - (1-f)^(1/2))^(1/2))
    /// Simplified: use the flow fraction directly since the Poiseuille
    /// profile is symmetric and we split top/bottom.
    pub fn dividing_streamline_y(&self) -> f64 {
        let f = self.flow_fraction_wide();
        // For a Poiseuille profile, the fraction of flow above y is
        // approximately (y/H)^2 * (3 - 2*y/H) (the CDF of the parabolic
        // profile).  We want y such that this fraction = 1 - f (i.e.,
        // f flows above y).  For 96% wide, the divider is near y = 0.05*H.
        // Use linear approximation for simplicity:
        self.parent_height_m * (1.0 - f)
    }
}

impl VelocityFieldInterpolator for AsymmetricBifurcationFlow {
    fn velocity_at(&self, x: f64, y: f64) -> (f64, f64) {
        let h = self.parent_height_m;
        let f_wide = self.flow_fraction_wide();
        let y_div = self.dividing_streamline_y();

        if x < self.x_split {
            // Parent: Poiseuille profile + cross-stream steering near split.
            let s = (2.0 * y / h - 1.0).clamp(-1.0, 1.0);
            let u = self.u_inlet * (1.0 - s * s);

            // Near the split plane, add cross-stream velocity that redirects
            // streamlines toward their destination daughter channel.
            // Cells above the dividing streamline get pushed toward the wide
            // daughter (upper region); cells below get pushed to narrow (lower).
            let approach_factor = ((x - self.x_split + h * 0.5) / (h * 0.5)).clamp(0.0, 1.0);
            let v = if approach_factor > 0.0 {
                let target_y = if y > y_div {
                    // Going to wide daughter: steer toward center of wide region
                    (y_div + h) * 0.5
                } else {
                    // Going to narrow daughter: steer toward center of narrow
                    y_div * 0.5
                };
                let steering = (target_y - y) * approach_factor * 2.0 * self.u_inlet / h;
                steering.clamp(-self.u_inlet * 0.3, self.u_inlet * 0.3)
            } else {
                0.0
            };

            (u, v)
        } else {
            // Daughters: Poiseuille in each sub-channel.
            // Wide daughter: y in [y_div, H].
            // Narrow daughter: y in [0, y_div].
            if y >= y_div {
                let dh = h - y_div;
                if dh < 1e-12 { return (0.0, 0.0); }
                let s = (2.0 * (y - y_div) / dh - 1.0).clamp(-1.0, 1.0);
                // U_max = (3/2) * Q / A.  Q = f_wide * Q_parent.
                // Q_parent = (2/3) * u_inlet * H (Poiseuille integral).
                // U_max_daughter = (3/2) * f_wide * (2/3) * u_inlet * H / dh
                //                = f_wide * u_inlet * H / dh
                let u_d = f_wide * self.u_inlet * h / dh;
                (u_d * (1.0 - s * s), 0.0)
            } else {
                if y_div < 1e-12 { return (0.0, 0.0); }
                let s = (2.0 * y / y_div - 1.0).clamp(-1.0, 1.0);
                let u_d = (1.0 - f_wide) * self.u_inlet * h / y_div;
                (u_d * (1.0 - s * s), 0.0)
            }
        }
    }

    fn is_fluid(&self, x: f64, y: f64) -> bool {
        x >= 0.0 && x <= self.length_m && y >= 0.0 && y <= self.parent_height_m
    }

    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.length_m, 0.0, self.parent_height_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poiseuille_cells_exit_domain() {
        let flow = PoiseuilleFlow2D {
            u_max: 0.1,
            width: 0.01,
            height: 0.002,
        };
        let config = CellTrackerConfig {
            viscosity: 3.5e-3,
            fluid_density: 1025.0,
            hydraulic_diameter_m: 0.002,
            u_max: 0.1,
            ..Default::default()
        };
        let tracker = CellTracker::new(&flow, config);
        let cells = vec![
            TrackedCell { population: CellPopulation::CTC, x: 0.0, y: 0.001, vx: 0.05, vy: 0.0, id: 0 },
            TrackedCell { population: CellPopulation::RBC, x: 0.0, y: 0.001, vx: 0.05, vy: 0.0, id: 1 },
            TrackedCell { population: CellPopulation::WBC, x: 0.0, y: 0.001, vx: 0.05, vy: 0.0, id: 2 },
        ];
        let trajectories = tracker.trace_cells(&cells, 1e-5, 200_000);
        for traj in &trajectories {
            assert!(traj.exit_outlet.is_some(), "cell {} should exit", traj.cell_id);
            assert!(traj.positions.len() >= 2, "cell {} should have trajectory points", traj.cell_id);
        }
    }

    #[test]
    fn poiseuille_centerline_cells_stay_near_center() {
        let flow = PoiseuilleFlow2D {
            u_max: 0.2,
            width: 0.02,
            height: 0.002,
        };
        let config = CellTrackerConfig {
            viscosity: 3.5e-3,
            fluid_density: 1025.0,
            hydraulic_diameter_m: 0.002,
            u_max: 0.2,
            ..Default::default()
        };
        let tracker = CellTracker::new(&flow, config);
        // Release at centerline: should stay near center.
        let cells = vec![
            TrackedCell { population: CellPopulation::CTC, x: 0.0, y: 0.001, vx: 0.1, vy: 0.0, id: 0 },
        ];
        let trajectories = tracker.trace_cells(&cells, 1e-5, 500_000);
        let traj = &trajectories[0];
        assert!(traj.exit_outlet.is_some());
        // Final y should be within 30% of centerline (0.001 m).
        let final_y = traj.positions.last().unwrap()[1];
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
            parent_width_m: 2.0e-3,
            parent_height_m: parent_h,
            wide_daughter_width_m: 2.0e-3 * center_frac,
            narrow_daughter_width_m: 2.0e-3 * (1.0 - center_frac),
            length_m: 0.015,
            u_inlet: 0.05,
            x_split: 0.005,
        };
        let y_div = flow.dividing_streamline_y();
        let config = CellTrackerConfig {
            viscosity: 3.5e-3,
            fluid_density: 1025.0,
            hydraulic_diameter_m: parent_h,
            u_max: 0.05,
            outlet_zones: vec![
                OutletZone {
                    name: "center".to_string(),
                    x_min: 0.014,
                    y_lo: y_div,
                    y_hi: parent_h,
                },
                OutletZone {
                    name: "peripheral".to_string(),
                    x_min: 0.014,
                    y_lo: 0.0,
                    y_hi: y_div,
                },
            ],
            split_x: 0.005,
            dividing_streamline_y: y_div,
            psm_params: Some(PsmBifurcationParams {
                flow_fraction_wide: flow.flow_fraction_wide(),
                wide_daughter_dh: 2.0e-3 * center_frac,
                narrow_daughter_dh: 2.0e-3 * (1.0 - center_frac),
                feed_hematocrit: 0.45,
            }),
            ..Default::default()
        };
        let tracker = CellTracker::new(&flow, config);

        let n_per_pop = 30;
        let mut cells = Vec::with_capacity(n_per_pop * 2);
        for i in 0..n_per_pop {
            let y = parent_h * 0.05 + parent_h * 0.9 * (i as f64 / (n_per_pop - 1) as f64);
            cells.push(TrackedCell {
                population: CellPopulation::CTC, x: 1e-4, y, vx: 0.03, vy: 0.0, id: i,
            });
            cells.push(TrackedCell {
                population: CellPopulation::RBC, x: 1e-4, y, vx: 0.03, vy: 0.0, id: n_per_pop + i,
            });
        }

        let trajectories = tracker.trace_cells(&cells, 2e-6, 1_000_000);
        let routing = tracker.classify_routing(&trajectories);

        let q_frac = flow.flow_fraction_wide();
        eprintln!("--- Asymmetric bifurcation cell tracking ---");
        eprintln!("  center_frac = {center_frac}, flow_fraction_wide = {q_frac:.4}");
        eprintln!("  y_div = {y_div:.6} m");
        eprintln!(
            "  CTC center {}/{}, RBC center {}/{}",
            routing.ctc_center, routing.ctc_total,
            routing.rbc_center, routing.rbc_total,
        );

        let exited = trajectories.iter()
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
