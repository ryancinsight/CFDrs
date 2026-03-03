//! DES length scale computations
//!
//! Contains the length scale computation methods for DES97, DDES, and IDDES variants.
//!
//! # DES97
//! L_DES = min(L_RANS, C_DES * Δ)
//!
//! # DDES (Spalart et al. 2006)
//! Uses shielding function: f_d = 1 - tanh[(8 r_d)^3]
//! where r_d = (ν_t + ν) / (κ² d_w² S)
//!
//! # IDDES (Shur et al. 2008)
//! Combines DDES shielding with WMLES capability via blending function f_B.

use super::{DESVariant, DetachedEddySimulation};
use nalgebra::DMatrix;

impl DetachedEddySimulation {
    /// Compute DES length scale
    pub(super) fn compute_des_length_scale(
        &self,
        velocity_u: &DMatrix<f64>,
        velocity_v: &DMatrix<f64>,
        dx: f64,
        dy: f64,
    ) -> DMatrix<f64> {
        let nx = velocity_u.nrows();
        let ny = velocity_v.ncols();
        let mut length_scale = DMatrix::zeros(nx, ny);

        // Compute strain rate magnitude for DDES and IDDES
        // Needed for shielding function r_d
        let strain_magnitude =
            if matches!(self.config.variant, DESVariant::DDES | DESVariant::IDDES) {
                self.compute_strain_rate_magnitude(velocity_u, velocity_v, dx, dy)
            } else {
                DMatrix::zeros(1, 1) // Dummy for DES97
            };

        for i in 0..nx {
            for j in 0..ny {
                // Grid spacing (max of dx, dy for simplicity)
                let delta = dx.max(dy);

                // RANS length scale is the wall distance in SA model
                let rans_length = self.wall_distance[(i, j)];

                // Get current modified viscosity state (nu_tilde)
                let nu_tilde_val = self.nu_tilde[(i, j)];

                // DES length scale based on variant
                let des_length = match self.config.variant {
                    DESVariant::DES97 => {
                        // Original DES: L_DES = min(L_RANS, C_DES * Δ)
                        rans_length.min(self.config.des_constant * delta)
                    }
                    DESVariant::DDES => {
                        // Delayed DES with shielding function
                        let strain_mag = strain_magnitude[(i, j)];
                        self.compute_ddes_length_scale(
                            rans_length,
                            delta,
                            nu_tilde_val,
                            strain_mag,
                            i,
                            j,
                        )
                    }
                    DESVariant::IDDES => {
                        // Improved DDES
                        let strain_mag = strain_magnitude[(i, j)];
                        self.compute_iddes_length_scale(
                            rans_length,
                            dx,
                            dy,
                            nu_tilde_val,
                            strain_mag,
                            i,
                            j,
                        )
                    }
                };

                length_scale[(i, j)] = des_length;
            }
        }

        length_scale
    }

    /// Compute DDES length scale with shielding function
    fn compute_ddes_length_scale(
        &self,
        rans_length: f64,
        delta: f64,
        nu_tilde_val: f64,
        strain_mag: f64,
        i: usize,
        j: usize,
    ) -> f64 {
        let d_w = self.wall_distance[(i, j)];

        match self.config.variant {
            DESVariant::DDES => {
                let l_rans = rans_length;
                let l_les = self.config.des_constant * delta;

                if l_rans < l_les {
                    // RANS mode - no shielding needed
                    l_rans
                } else {
                    // LES mode with proper DDES shielding function
                    // Spalart et al. (2006): fd = 1 - tanh[(8 * r_d)^3]
                    // r_d = (nu_t + nu) / (kappa^2 * d_w^2 * max(S, 1e-10))

                    let kappa = 0.41;

                    // Compute turbulent viscosity from nu_tilde using SA relation
                    let nu_t = self
                        .spalart_allmaras
                        .eddy_viscosity(nu_tilde_val, self.config.rans_viscosity);
                    let nu = self.config.rans_viscosity;

                    // Strain rate magnitude
                    let s = strain_mag.max(1e-10);

                    // r_d calculation
                    let denom = kappa * kappa * d_w * d_w * s;
                    let r_d = (nu_t + nu) / denom.max(1e-20);

                    // Shielding function
                    let f_d = 1.0 - (8.0 * r_d).powi(3).tanh();

                    // Blend RANS and LES length scales
                    // f_d is 0 in RANS region (r_d large), 1 in LES region (r_d small)
                    l_rans * (1.0 - f_d) + l_les * f_d
                }
            }
            _ => rans_length.min(self.config.des_constant * delta),
        }
    }

    /// Compute IDDES length scale
    fn compute_iddes_length_scale(
        &self,
        rans_length: f64,
        dx: f64,
        dy: f64,
        nu_tilde_val: f64,
        strain_mag: f64,
        i: usize,
        j: usize,
    ) -> f64 {
        // IDDES formulation (Shur et al. 2008, Gritskevich et al. 2012)
        // Combines DDES with WMLES capabilities.

        // Constants
        let c_w = 0.15;
        let c_des = self.config.des_constant;
        let kappa = 0.41;

        // Grid scales
        let h_max = dx.max(dy);
        let h_wn = dx.min(dy); // Approximation for wall-normal spacing

        // Wall distance
        let d_w = self.wall_distance[(i, j)];

        // 1. Definition of Delta_IDDES
        // Delta_IDDES = min(max(Cw * dw, Cw * h_max, h_wn), h_max)
        let delta_iddes = (c_w * d_w).max(c_w * h_max).max(h_wn).min(h_max);

        // 2. LES length scale
        let l_les = c_des * delta_iddes;

        // 3. Shielding function f_dt (similar to DDES) and blending f_B
        // r_dt = nu_t / (k^2 * dw^2 * S)

        // Compute actual turbulent viscosity
        let nu_t = self
            .spalart_allmaras
            .eddy_viscosity(nu_tilde_val, self.config.rans_viscosity);

        let s = strain_mag.max(1e-10);

        // Use actual nu_t instead of mixing length approximation
        let r_dt = nu_t / (kappa * kappa * d_w * d_w * s).max(1e-20);

        // f_dt = 1 - tanh[(8 * r_dt)^3] (DDES-like shielding)
        // Using 8.0 to match standard DDES formulation.
        let f_dt = 1.0 - (8.0 * r_dt).powi(3).tanh();

        // f_B = tanh[(16 * r_dt)^3] (Blending function)
        // Using 16.0 to ensure f_B rises before f_dt drops, maintaining RANS shielding
        // in the transition region.
        let f_b = (16.0 * r_dt).powi(3).tanh();

        // Final shielding function
        // tilde_f_d = max(1 - f_dt, f_B)
        let tilde_f_d = (1.0 - f_dt).max(f_b);

        // 4. IDDES length scale
        // l_IDDES = tilde_f_d * l_RANS + (1 - tilde_f_d) * l_LES
        // Note: Omitting elevating function f_e for standard hybrid behavior.

        tilde_f_d * rans_length + (1.0 - tilde_f_d) * l_les
    }

    /// Compute strain rate magnitude (shared with Smagorinsky)
    pub(super) fn compute_strain_rate_magnitude(
        &self,
        velocity_u: &DMatrix<f64>,
        velocity_v: &DMatrix<f64>,
        dx: f64,
        dy: f64,
    ) -> DMatrix<f64> {
        let nx = velocity_u.nrows();
        let ny = velocity_u.ncols();
        let mut strain_magnitude = DMatrix::zeros(nx, ny);

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                // Velocity gradients
                let du_dx = (velocity_u[(i + 1, j)] - velocity_u[(i - 1, j)]) / (2.0 * dx);
                let du_dy = (velocity_u[(i, j + 1)] - velocity_u[(i, j - 1)]) / (2.0 * dy);
                let dv_dx = (velocity_v[(i + 1, j)] - velocity_v[(i - 1, j)]) / (2.0 * dx);
                let dv_dy = (velocity_v[(i, j + 1)] - velocity_v[(i, j - 1)]) / (2.0 * dy);

                // Strain rate tensor magnitude
                let s11 = du_dx;
                let s22 = dv_dy;
                let s12 = 0.5 * (du_dy + dv_dx);

                strain_magnitude[(i, j)] =
                    (2.0 * s11 * s11 + 2.0 * s22 * s22 + 4.0 * s12 * s12).sqrt();
            }
        }

        strain_magnitude
    }
}
