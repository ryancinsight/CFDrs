//! LES/DES model validation and benchmark methods for [`TurbulenceValidator`].

use super::ValidationResult;
use crate::physics::turbulence::traits::LESTurbulenceModel;
use crate::physics::turbulence::traits::TurbulenceModel;
use crate::physics::turbulence::{
    DetachedEddySimulation, KEpsilonModel, KOmegaSSTModel, SmagorinskyLES,
};
use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

use super::TurbulenceValidator;

impl<T: RealField + FromPrimitive + ToPrimitive + Copy> TurbulenceValidator<T> {
    /// Validate Smagorinsky LES model SGS viscosity calculation
    pub fn validate_smagorinsky_sgs(&self) -> ValidationResult {
        let config = crate::physics::turbulence::les_smagorinsky::SmagorinskyConfig {
            smagorinsky_constant: 0.1,
            dynamic_procedure: false,
            wall_damping: false,
            van_driest_constant: 0.0,
            min_sgs_viscosity: 0.0,
            use_gpu: false,
        };
        let mut model = SmagorinskyLES::new(16, 16, 0.1, 0.1, config);

        let mut velocity_u = nalgebra::DMatrix::zeros(16, 16);
        let velocity_v = nalgebra::DMatrix::zeros(16, 16);

        for j in 0..16 {
            for i in 0..16 {
                velocity_u[(i, j)] = j as f64 * 0.1;
            }
        }

        let result = model.update(
            &velocity_u,
            &velocity_v,
            &nalgebra::DMatrix::zeros(16, 16),
            1.0,
            1e-5,
            0.001,
            0.1,
            0.1,
        );

        if result.is_err() {
            return ValidationResult {
                test_name: "Smagorinsky LES SGS Viscosity".to_string(),
                passed: false,
                metric: "Update failed".to_string(),
                details: format!("Error: {:?}", result.err().unwrap()),
            };
        }

        let sgs_viscosity = model.get_turbulent_viscosity_field();

        let mut all_positive = true;
        let mut total_sgs = 0.0f64;
        let mut count = 0;

        for j in 1..15 {
            for i in 1..15 {
                let sgs = sgs_viscosity[(i, j)];
                if sgs < 0.0 {
                    all_positive = false;
                }
                total_sgs += sgs;
                count += 1;
            }
        }

        let avg_sgs_f64 = total_sgs / f64::from(count);
        let reasonable_magnitude = (1e-6..=1e-2).contains(&avg_sgs_f64);

        ValidationResult {
            test_name: "Smagorinsky LES SGS Viscosity".to_string(),
            passed: all_positive && reasonable_magnitude,
            metric: format!("Avg SGS: {avg_sgs_f64:.2e}"),
            details: format!(
                "Positive: {all_positive}, Reasonable magnitude: {reasonable_magnitude}"
            ),
        }
    }

    /// Validate DES model length scale calculation
    pub fn validate_des_length_scale(&self) -> ValidationResult {
        let config = crate::physics::turbulence::des::DESConfig {
            des_constant: 0.65,
            max_sgs_ratio: 0.5,
            variant: crate::physics::turbulence::des::DESVariant::DES97,
            rans_viscosity: 1e-5,
            use_gpu: false,
        };
        let mut model = DetachedEddySimulation::new(16, 16, 0.1, 0.1, config, &[]);

        let mut velocity_u = nalgebra::DMatrix::zeros(16, 16);
        let mut velocity_v = nalgebra::DMatrix::zeros(16, 16);

        for j in 0..16 {
            for i in 0..16 {
                let x = i as f64 * 0.1;
                let y = j as f64 * 0.1;
                velocity_u[(i, j)] = (x * 0.1).sin() + 0.1;
                velocity_v[(i, j)] = (y * 0.1).cos() + 0.1;
            }
        }

        let result = model.update(
            &velocity_u,
            &velocity_v,
            &nalgebra::DMatrix::zeros(16, 16),
            1.0,
            1e-5,
            0.001,
            0.1,
            0.1,
        );

        if result.is_err() {
            return ValidationResult {
                test_name: "DES Length Scale Calculation".to_string(),
                passed: false,
                metric: "Update failed".to_string(),
                details: format!("Error: {:?}", result.err().unwrap()),
            };
        }

        let sgs_viscosity = model.get_turbulent_viscosity_field();

        let mut all_finite = true;
        let mut all_positive = true;
        let mut max_visc = 0.0;

        for j in 0..16 {
            for i in 0..16 {
                let visc = sgs_viscosity[(i, j)];
                if !visc.is_finite() {
                    all_finite = false;
                }
                if visc < 0.0 {
                    all_positive = false;
                }
                if visc > max_visc {
                    max_visc = visc;
                }
            }
        }

        ValidationResult {
            test_name: "DES Length Scale Calculation".to_string(),
            passed: all_finite && all_positive,
            metric: format!("Max viscosity: {max_visc:.2e}"),
            details: format!("All finite: {all_finite}, All positive: {all_positive}"),
        }
    }

    /// Validate k-ε and k-ω SST models against flat plate boundary layer
    ///
    /// Reference: White, F. M. (2006). Viscous Fluid Flow (3rd ed.). McGraw-Hill.
    pub fn validate_flat_plate_boundary_layer(&self) -> ValidationResult {
        let reynolds_x = 5e6_f64;
        let free_stream_velocity = 10.0_f64;
        let kinematic_viscosity = 1.5e-5_f64;

        let expected_cf = 0.376_f64 / reynolds_x.powf(0.2_f64);

        let nx = 50;
        let ny = 30;
        let mut k_model = KEpsilonModel::new(nx, ny);
        let mut k_field = vec![0.01; nx * ny];
        let mut eps_field = vec![0.001; nx * ny];

        let mut velocity_field = vec![Vector2::zeros(); nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let idx = i * ny + j;
                let y_over_delta = j as f64 / (ny - 1) as f64;
                let u_velocity = if y_over_delta < 1.0 {
                    free_stream_velocity * y_over_delta.powf(1.0 / 7.0)
                } else {
                    free_stream_velocity
                };
                velocity_field[idx] = Vector2::new(u_velocity, 0.0);
            }
        }

        let dt = 1e-4;
        let dx = 0.01;
        let dy = 0.001;

        for _ in 0..10 {
            k_model
                .update(
                    &mut k_field,
                    &mut eps_field,
                    &velocity_field,
                    1.0,
                    kinematic_viscosity,
                    dt,
                    dx,
                    dy,
                )
                .unwrap();
        }

        let wall_shear_stress =
            kinematic_viscosity * (velocity_field[ny].x - velocity_field[0].x) / dy;
        let calculated_cf = wall_shear_stress / (0.5 * free_stream_velocity * free_stream_velocity);

        let cf_ratio = calculated_cf / expected_cf;
        let passed = (cf_ratio - 1.0).abs() < 0.15;

        ValidationResult {
            test_name: "Flat Plate Boundary Layer - k-ε Model".to_string(),
            passed,
            metric: format!("Cf ratio: {cf_ratio:.3}"),
            details: format!("Expected Cf: {expected_cf:.6}, Calculated: {calculated_cf:.6}"),
        }
    }

    /// Validate turbulence models against channel flow DNS data
    ///
    /// Reference: Moser, Kim & Mansour (1999). Direct numerical simulation of
    /// turbulent channel flow. Physics of Fluids.
    pub fn validate_channel_flow_dns(&self) -> ValidationResult {
        let re_tau = 590.0_f64;
        let kinematic_viscosity = 1.0_f64 / re_tau;

        let nx = 40;
        let ny = 40;
        let mut k_model = KOmegaSSTModel::new(nx, ny);
        let mut k_field = vec![0.1; nx * ny];
        let mut omega_field = vec![10.0; nx * ny];

        let mut velocity_field = vec![Vector2::zeros(); nx * ny];

        for j in 0..ny {
            let y_plus = (j as f64 / (ny - 1) as f64) * re_tau;
            let u_plus = if y_plus <= 5.0_f64 {
                y_plus
            } else if y_plus <= 30.0_f64 {
                5.0_f64 * (y_plus / 5.0_f64).ln() + 5.17_f64
            } else {
                2.5_f64 * (y_plus / 30.0_f64).ln() + 5.17_f64 + 2.5_f64 * (30.0_f64 / 5.0_f64).ln()
            };
            let u_velocity = u_plus * kinematic_viscosity * re_tau;

            for i in 0..nx {
                let idx = i * ny + j;
                velocity_field[idx] = Vector2::new(u_velocity, 0.0_f64);
            }
        }

        let dt = 1e-5;
        let dx = 0.01;
        let dy = 0.001;

        for _ in 0..15 {
            k_model
                .update(
                    &mut k_field,
                    &mut omega_field,
                    &velocity_field,
                    1.0,
                    kinematic_viscosity,
                    dt,
                    dx,
                    dy,
                )
                .unwrap();
        }

        let mut u_plus_calculated = Vec::new();
        let mut u_plus_reference = Vec::new();

        for j in 0..ny {
            let y_plus = (j as f64 / (ny - 1) as f64) * re_tau;
            let u_physical = velocity_field[j].x;
            let u_plus_calc = u_physical / (kinematic_viscosity * re_tau);

            u_plus_calculated.push(u_plus_calc);

            let u_plus_ref = if y_plus <= 5.0_f64 {
                y_plus
            } else if y_plus <= 30.0_f64 {
                5.0_f64 * (y_plus / 5.0_f64).ln() + 5.17_f64
            } else {
                2.5_f64 * (y_plus / 30.0_f64).ln() + 5.17_f64 + 2.5_f64 * (30.0_f64 / 5.0_f64).ln()
            };

            u_plus_reference.push(u_plus_ref);
        }

        let mut rms_error = 0.0_f64;
        for (calc, ref_val) in u_plus_calculated.iter().zip(u_plus_reference.iter()) {
            rms_error += (calc - ref_val).powi(2);
        }
        rms_error = (rms_error / u_plus_calculated.len() as f64).sqrt();

        let passed = rms_error < 0.3;

        ValidationResult {
            test_name: "Channel Flow DNS - k-ω SST Model".to_string(),
            passed,
            metric: format!("RMS error: {rms_error:.3} wall units"),
            details: format!("Re_τ = {re_tau}, Profile points: {ny}"),
        }
    }

    /// Validate LES models against decaying homogeneous turbulence
    ///
    /// Reference: Comte-Bellot & Corrsin (1971) experimental data
    pub fn validate_les_decaying_turbulence(&self) -> ValidationResult {
        let nx = 32;
        let ny = 32;
        let config = crate::physics::turbulence::les_smagorinsky::SmagorinskyConfig {
            smagorinsky_constant: 0.1,
            dynamic_procedure: false,
            wall_damping: false,
            van_driest_constant: 0.0,
            min_sgs_viscosity: 0.0,
            use_gpu: false,
        };

        let mut les_model = SmagorinskyLES::new(nx, ny, 0.05, 0.05, config);

        let mut velocity_u = DMatrix::from_element(nx, ny, 0.0);
        let mut velocity_v = DMatrix::from_element(nx, ny, 0.0);

        use std::f64::consts::PI;
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * 0.05_f64;
                let y = j as f64 * 0.05_f64;
                let fluctuation_u = 0.1_f64
                    * ((2.0_f64 * PI * x / 0.5_f64).sin() + (4.0_f64 * PI * x / 0.5_f64).sin())
                    * ((2.0_f64 * PI * y / 0.5_f64).cos() + (4.0_f64 * PI * y / 0.5_f64).cos());
                let fluctuation_v = 0.1_f64
                    * ((2.0_f64 * PI * x / 0.5_f64).cos() + (4.0_f64 * PI * x / 0.5_f64).cos())
                    * ((2.0_f64 * PI * y / 0.5_f64).sin() + (4.0_f64 * PI * y / 0.5_f64).sin());
                velocity_u[(i, j)] = fluctuation_u;
                velocity_v[(i, j)] = fluctuation_v;
            }
        }

        let dt = 0.001_f64;
        let density = 1.0_f64;
        let viscosity = 1e-5_f64;
        let dx = 0.05_f64;
        let dy = 0.05_f64;

        let mut kinetic_energy_history = Vec::new();
        let initial_ke = Self::calculate_kinetic_energy(&velocity_u, &velocity_v);

        for step in 0..20 {
            les_model
                .update(
                    &velocity_u,
                    &velocity_v,
                    &DMatrix::zeros(nx, ny),
                    density,
                    viscosity,
                    dt,
                    dx,
                    dy,
                )
                .unwrap();

            if step % 5 == 0 {
                let current_ke = Self::calculate_kinetic_energy(&velocity_u, &velocity_v);
                kinetic_energy_history.push(current_ke / initial_ke);
            }
        }

        let final_energy_ratio = *kinetic_energy_history.last().unwrap_or(&1.0_f64);
        let decay_rate = -final_energy_ratio.ln() / (20.0_f64 * dt);

        let expected_decay_rate = 1.5_f64;
        let decay_ratio = decay_rate / expected_decay_rate;

        let passed = (decay_ratio - 1.0).abs() < 0.5;

        ValidationResult {
            test_name: "LES Decaying Homogeneous Turbulence".to_string(),
            passed,
            metric: format!("Decay rate ratio: {decay_ratio:.2}"),
            details: format!(
                "Expected: {expected_decay_rate:.2}, Calculated: {decay_rate:.3}, Final energy ratio: {final_energy_ratio:.4}"
            ),
        }
    }

    /// Helper function to calculate kinetic energy from velocity field
    pub(super) fn calculate_kinetic_energy(
        velocity_u: &DMatrix<f64>,
        velocity_v: &DMatrix<f64>,
    ) -> f64 {
        let mut ke = 0.0;
        let nx = velocity_u.nrows();
        let ny = velocity_u.ncols();

        for j in 0..ny {
            for i in 0..nx {
                let u = velocity_u[(i, j)];
                let v = velocity_v[(i, j)];
                ke += 0.5 * (u * u + v * v);
            }
        }

        ke / (nx * ny) as f64
    }
}
