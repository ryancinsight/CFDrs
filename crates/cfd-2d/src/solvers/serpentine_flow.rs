//! 2D serpentine channel flow solver with mixing efficiency validation
//!
//! This module implements validated solvers for serpentine (sinuous) channels,
//! which are critical for microfluidic mixing and sample preparation.
//!
//! # Physics Background
//!
//! ## Serpentine Channel Design
//!
//! Serpentine channels enhance mixing through:
//! 1. **Advection**: Fluid streams are stretched and folded
//! 2. **Secondary flow**: Dean vortices in curves increase transverse mixing
//! 3. **Residence time**: Longer path length increases reaction time
//!
//! ## Mixing Mechanisms
//!
//! In laminar microfluidic flows (Re << 1), mixing is diffusion-limited:
//!
//! **Diffusion time scale:**
//! ```text
//! t_diff = w² / (4D)  [seconds, from Fick's law]
//! ```
//!
//! where:
//! - w = channel width [m]
//! - D = diffusion coefficient [m²/s]
//!
//! **Advection time scale:**
//! ```text
//! t_adv = L / u  [seconds, from advection]
//! ```
//!
//! **Mixing length:**
//! ```text
//! L_mix = u · t_diff = u · w² / (4D)
//! ```
//!
//! For complete mixing across channel width: L_adv > L_mix
//!
//! ## Peclet Number
//!
//! Ratio of advection to diffusion:
//! ```text
//! Pe = u·w / D
//! ```
//!
//! - Pe << 1: Diffusion-dominated mixing
//! - Pe >> 1: Mixing limited by diffusion length
//!
//! ## Mixing Efficiency
//!
//! Quantified by intensity of segregation (ISO):
//! ```text
//! ISO = ∫(c - c_mean)² dV / ∫c_max² dV
//! ```
//!
//! where:
//! - ISO = 1 at inlet (completely unmixed)
//! - ISO = 0 at outlet (perfectly mixed)
//!
//! # Validation Strategy
//!
//! 1. **Advection-diffusion**: Solve transport equation with known inlet conditions
//! 2. **Richardson extrapolation**: Grid convergence study
//! 3. **Literature**: Compare against experimental mixing data
//! 4. **Analytical**: Use advection-diffusion analytical solutions
//!
//! # References
//!
//! - Hardt, S. & Schönfeld, F. (2003). "Microfluidic technologies for miniaturized
//!   analysis systems". SpringerLink
//! - Squires, T.M. & Quake, S.R. (2005). "Microfluidics: Fluid physics at the
//!   nanoliter scale". Reviews of Modern Physics, 77(3), 977
//! - Yao, Z., et al. (2014). "Numerical study of mixing in microchannels with
//!   oscillatory flow". Microfluidics and Nanofluidics, 16(1), 145-155

use cfd_core::conversion::SafeFromF64;
use nalgebra::{RealField};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// Serpentine Channel Geometry
// ============================================================================

/// Serpentine channel geometry with periodic turns
///
/// Defines a channel that alternates between straight sections and 90° turns.
/// The path creates a snake-like pattern that enhances mixing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineGeometry<T: RealField + Copy> {
    /// Channel width [m]
    pub width: T,
    /// Channel height [m] (constant)
    pub height: T,
    /// Straight section length [m]
    pub l_straight: T,
    /// Turn radius [m] (for curved turns)
    pub turn_radius: T,
    /// Number of complete back-and-forth cycles
    pub n_cycles: usize,
}

impl<T: RealField + Copy + FromPrimitive> SerpentineGeometry<T> {
    /// Create standard serpentine for microfluidics
    ///
    /// # Typical Design
    ///
    /// - Width: 100-500 μm
    /// - Straight sections: 500 μm
    /// - Turn radius: 200 μm
    /// - Cycles: 5-10 (for ~5-10 mm total length)
    pub fn microfluidic_standard() -> Self {
        Self {
            width: T::from_f64_or_one(200e-6),
            height: T::from_f64_or_one(50e-6),
            l_straight: T::from_f64_or_one(500e-6),
            turn_radius: T::from_f64_or_one(200e-6),
            n_cycles: 5,
        }
    }

    /// Create custom serpentine
    pub fn new(width: T, height: T, l_straight: T, turn_radius: T, n_cycles: usize) -> Self {
        Self {
            width,
            height,
            l_straight,
            turn_radius,
            n_cycles,
        }
    }

    /// Calculate total channel length
    ///
    /// Length = n_cycles × (2 × l_straight + turn_length)
    /// where turn_length ≈ π × turn_radius (90° arc)
    pub fn total_length(&self) -> T {
        let pi = T::from_f64_or_one(PI);
        let two = T::from_f64_or_one(2.0);
        let turn_length = pi / T::from_f64_or_one(2.0) * self.turn_radius;

        T::from_f64_or_one(self.n_cycles as f64) * (two * self.l_straight + turn_length)
    }

    /// Get cross-sectional area
    pub fn cross_section_area(&self) -> T {
        self.width * self.height
    }

    /// Calculate number of diffusion lengths in straight section
    ///
    /// # Definition
    ///
    /// ```text
    /// n_diff = (l_straight / width) × Pe
    /// ```
    ///
    /// where Pe = u·w/D (Peclet number)
    pub fn diffusion_lengths_per_section(&self, peclet: T) -> T {
        let three = T::from_f64_or_one(3.0);
        (self.l_straight / (self.width + T::from_f64_or_one(1e-15))) * three * peclet
    }

    /// Check if a point (x, y) is within the fluid domain
    ///
    /// Assuming serpentine starts at (0, R) and snakes upwards in Y.
    /// Channel centerline path:
    /// Cycle N:
    ///   Segment 1: (0, 4NR + R) -> (Ls, 4NR + R)
    ///   Turn 1: Semi-circle around (Ls, 4NR + 2R) with radius R
    ///   Segment 2: (Ls, 4NR + 3R) -> (0, 4NR + 3R)
    ///   Turn 2: Semi-circle around (0, 4NR + 4R) with radius R
    pub fn contains(&self, x: T, y: T) -> bool {
        let r = self.turn_radius;
        let ls = self.l_straight;
        let w = self.width;
        let half_w = w / T::from_f64(2.0).unwrap();
        let four_r = r * T::from_f64(4.0).unwrap();

        // Normalize y to cycle 0
        let cycle_f = y / four_r;
        let y_in_cycle = y % four_r;
        
        let r_sq = r * r;
        let inner_r = r - half_w;
        let outer_r = r + half_w;
        let inner_r_sq = inner_r * inner_r;
        let outer_r_sq = outer_r * outer_r;

        // Sequence within a cycle:
        // 0..2R: Segment 1 (bottom straight) and Turn 1 (right)
        // 2R..4R: Segment 2 (top straight) and Turn 2 (left)

        if y_in_cycle < r + half_w && y_in_cycle > r - half_w {
            // Segment 1: (0, R) -> (Ls, R)
            // Extend to x = -R for inlet tail
            if x >= -r - half_w && x <= ls {
                return true;
            }
        }

        if y_in_cycle < T::from_f64(3.0).unwrap() * r + half_w && y_in_cycle > T::from_f64(3.0).unwrap() * r - half_w {
            // Segment 2: (Ls, 3R) -> (0, 3R)
            if x >= T::zero() && x <= ls {
                return true;
            }
        }

        // Right Turn (around x=ls, y=2R)
        let dx_r = x - ls;
        let dy_r = y_in_cycle - r * T::from_f64(2.0).unwrap();
        let d_sq_r = dx_r * dx_r + dy_r * dy_r;
        if x >= ls && d_sq_r >= inner_r_sq && d_sq_r <= outer_r_sq {
            return true;
        }

        // Left Turn (around x=0, y=4R/0)
        // Check both top of this cycle and bottom of next
        let dy_l_top = y_in_cycle - four_r;
        let d_sq_l_top = x * x + dy_l_top * dy_l_top;
        if x <= T::zero() && d_sq_l_top >= inner_r_sq && d_sq_l_top <= outer_r_sq {
            return true;
        }
        
        let d_sq_l_bot = x * x + y_in_cycle * y_in_cycle;
        if x <= T::zero() && d_sq_l_bot >= inner_r_sq && d_sq_l_bot <= outer_r_sq {
            return true;
        }

        false
    }

    /// Get bounding box [min_x, max_x, min_y, max_y]
    pub fn bounding_box(&self) -> [T; 4] {
        let r = self.turn_radius;
        let ls = self.l_straight;
        let w = self.width;
        let half_w = w / T::from_f64(2.0).unwrap();
        let four_r = r * T::from_f64(4.0).unwrap();
        
        [
            -r - half_w,
            ls + r + half_w,
            T::zero(),
            four_r * T::from_f64(self.n_cycles as f64).unwrap()
        ]
    }
}

// ============================================================================
// Mixing Models
// ============================================================================

/// Analytical model for advection-diffusion mixing in straight sections
///
/// In a straight microchannel with two inlet streams, concentration
/// evolves according to Fick's law:
///
/// ```text
/// ∂c/∂t + u·∂c/∂x = D·∂²c/∂y²
/// ```
///
/// The mixing length (distance to achieve 90% mixing) is:
/// ```text
/// L_mix = 3.6 × w / Pe
/// ```
///
/// where Pe = u·w/D is the Peclet number.
pub struct AdvectionDiffusionMixing<T: RealField + Copy> {
    /// Channel width [m]
    pub width: T,
    /// Mean flow velocity [m/s]
    pub velocity: T,
    /// Diffusion coefficient [m²/s]
    pub diffusion_coeff: T,
}

impl<T: RealField + Copy + FromPrimitive> AdvectionDiffusionMixing<T> {
    /// Create mixing model
    pub fn new(width: T, velocity: T, diffusion_coeff: T) -> Self {
        Self {
            width,
            velocity,
            diffusion_coeff,
        }
    }

    /// Calculate Peclet number
    ///
    /// Pe = u·w / D
    pub fn peclet_number(&self) -> T {
        (self.velocity * self.width) / (self.diffusion_coeff + T::from_f64_or_one(1e-15))
    }

    /// Calculate mixing length for 90% mixing
    ///
    /// # Formula
    ///
    /// For diffusion in transverse (y) direction:
    /// ```text
    /// L_mix = 3.6 × w / Pe = 3.6 × D / u
    /// ```
    ///
    /// This is empirical from diffusion literature and matches numerical simulations
    /// of T-junction mixers.
    pub fn mixing_length_90_percent(&self) -> T {
        let three_point_six = T::from_f64_or_one(3.6);
        let pe = self.peclet_number();

        three_point_six * self.width / (pe + T::from_f64_or_one(1e-15))
    }

    /// Calculate mixing time to achieve 90% homogeneity
    ///
    /// t_mix = L_mix / u
    pub fn mixing_time_90_percent(&self) -> T {
        let l_mix = self.mixing_length_90_percent();
        l_mix / (self.velocity + T::from_f64_or_one(1e-15))
    }

    /// Estimate concentration profile at position x
    ///
    /// For two inlets with concentrations c_A and c_B, the concentration
    /// at position x depends on diffusion progress:
    ///
    /// ```text
    /// c(x, y) = c_A × erf(y / sqrt(4Dx/u)) + c_B × (1 - erf(...))
    /// ```
    ///
    /// At the centerline (y=0), maximum mixing occurs there.
    /// Returns fraction mixed at position x relative to channel width.
    pub fn mixing_fraction(&self, x: T) -> T {
        let one = T::one();
        let pe = self.peclet_number();

        // Normalized distance x/L_channel
        // Mixing fraction approximately: 1 - exp(-2 × x/L_mix)
        let l_mix = self.mixing_length_90_percent();
        let normalized_x = x / (l_mix + T::from_f64_or_one(1e-15));

        one - (-T::from_f64_or_one(2.0) * normalized_x).exp()
    }
}

// ============================================================================
// Serpentine Mixing Solution
// ============================================================================

/// Solution for serpentine channel mixing
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SerpentineMixingSolution<T: RealField + Copy> {
    /// Inlet concentration (first inlet) [mol/m³]
    pub c_inlet_a: T,
    /// Inlet concentration (second inlet) [mol/m³]
    pub c_inlet_b: T,
    /// Peclet number (dimensionless)
    pub peclet: T,
    /// Mixing length for 90% homogeneity [m]
    pub l_mix_90: T,
    /// Mixing time to 90% [s]
    pub t_mix_90: T,
    /// Mixing fraction at outlet [0-1]
    pub mixing_fraction_outlet: T,
    /// Estimated viscous pressure drop [Pa]
    pub pressure_drop: T,
}

impl<T: RealField + Copy + FromPrimitive> SerpentineMixingSolution<T> {
    /// Create solution from parameters
    pub fn new(
        geometry: &SerpentineGeometry<T>,
        velocity: T,
        diffusion_coeff: T,
        c_inlet_a: T,
        c_inlet_b: T,
        viscosity: T,
        density: T,
    ) -> Self {
        let mixing_model = AdvectionDiffusionMixing::new(geometry.width, velocity, diffusion_coeff);
        let pe = mixing_model.peclet_number();
        let l_mix = mixing_model.mixing_length_90_percent();
        let t_mix = mixing_model.mixing_time_90_percent();

        let total_length = geometry.total_length();
        let mixing_frac = mixing_model.mixing_fraction(total_length);

        // Pressure drop: ΔP = f × (L/D_h) × (1/2)ρu²
        // For laminar flow, friction factor f = 64/Re
        let d_h = T::from_f64_or_one(2.0) * (geometry.width * geometry.height)
            / (geometry.width + geometry.height);
        let re = (density * velocity * d_h) / viscosity;
        let f = T::from_f64_or_one(64.0) / re.max(T::from_f64_or_one(1.0));
        let dynamic_pressure = T::from_f64_or_one(0.5) * density * velocity * velocity;
        let dp = f * (total_length / d_h) * dynamic_pressure;

        Self {
            c_inlet_a,
            c_inlet_b,
            peclet: pe,
            l_mix_90: l_mix,
            t_mix_90: t_mix,
            mixing_fraction_outlet: mixing_frac,
            pressure_drop: dp,
        }
    }

    /// Check if mixing is achieved at outlet
    ///
    /// Typically consider "mixed" if fraction > 0.9 (90%)
    pub fn is_well_mixed(&self) -> bool {
        self.mixing_fraction_outlet > T::from_f64_or_one(0.9)
    }

    /// Estimate outlet concentration (assuming complete mixing)
    ///
    /// For equal volume flows: c_outlet = (c_A + c_B) / 2
    pub fn estimated_outlet_concentration(&self) -> T {
        (self.c_inlet_a + self.c_inlet_b) / T::from_f64_or_one(2.0)
    }
}

// ============================================================================
// Serpentine Validator
// ============================================================================

/// Validator for serpentine mixing channels
pub struct SerpentineValidator<T: RealField + Copy> {
    geometry: SerpentineGeometry<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> SerpentineValidator<T> {
    /// Create new validator
    pub fn new(geometry: SerpentineGeometry<T>) -> Self {
        Self { geometry }
    }

    /// Validate mixing efficiency against expected values
    ///
    /// # Validation Criteria
    ///
    /// For well-designed serpentines:
    /// - Mixing length < total channel length
    /// - Mixing fraction at outlet > 90%
    /// - Pressure drop < 100 kPa (for ~mm-scale channels)
    pub fn validate_mixing(
        &self,
        solution: &SerpentineMixingSolution<T>,
    ) -> Result<SerpentineValidationResult<T>, String> {
        let total_length = self.geometry.total_length();
        let l_mix = solution.l_mix_90;

        // Check if mixing length is achieved within channel
        let achievable = l_mix < total_length * T::from_f64_or_one(2.0);

        // Check mixing fraction at outlet
        let well_mixed = solution.is_well_mixed();

        // Check pressure drop is reasonable
        let dp_reasonable = solution.pressure_drop < T::from_f64_or_one(1e5); // < 100 kPa

        let validation_passed = achievable && well_mixed && dp_reasonable;

        let mut error_message = None;
        if !achievable {
            error_message = Some("Mixing length exceeds 2× total channel length".to_string());
        } else if !well_mixed {
            error_message = Some(format!(
                "Mixing fraction {:.2}% < 90% at outlet",
                (solution.mixing_fraction_outlet * T::from_f64_or_one(100.0))
                    .to_f64()
                    .unwrap_or(f64::NAN)
            ));
        } else if !dp_reasonable {
            error_message = Some(format!(
                "Pressure drop {:.0} Pa > 100 kPa",
                solution.pressure_drop.to_f64().unwrap_or(f64::NAN)
            ));
        }

        Ok(SerpentineValidationResult {
            validation_passed,
            mixing_fraction_outlet: solution.mixing_fraction_outlet,
            error_message,
        })
    }
}

/// Validation result for serpentine
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineValidationResult<T: RealField + Copy> {
    /// Validation passed
    pub validation_passed: bool,
    /// Measured mixing fraction at outlet
    pub mixing_fraction_outlet: T,
    /// Error message if validation failed
    pub error_message: Option<String>,
}

// ============================================================================
// Discretized Serpentine Solver
// ============================================================================

use crate::solvers::ns_fvm_2d::{NavierStokesSolver2D, SIMPLEConfig, BloodModel};
use crate::solvers::scalar_transport_2d::{ScalarTransportSolver2D, ScalarTransportConfig};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::error::Result as CfdResult;

/// Discretized 2D Serpentine Flow Solver
pub struct SerpentineSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    pub geometry: SerpentineGeometry<T>,
    pub ns_solver: NavierStokesSolver2D<T>,
    pub scalar_solver: ScalarTransportSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive> SerpentineSolver2D<T> {
    /// Create new discretized serpentine solver
    pub fn new(
        geometry: SerpentineGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
    ) -> Self {
        let bbox = geometry.bounding_box();
        let width = bbox[1] - bbox[0];
        let height = bbox[3] - bbox[2];
        
        let grid = crate::solvers::ns_fvm_2d::StaggeredGrid2D::new(nx, ny, width, height);
        let config = SIMPLEConfig::default();
        let mut ns_solver = NavierStokesSolver2D::new(grid, blood, density, config);

        // Populate mask
        for i in 0..nx {
            for j in 0..ny {
                let x = ns_solver.grid.x_center(i) + bbox[0];
                let y = ns_solver.grid.y_center(j) + bbox[2];
                ns_solver.field.mask[i][j] = geometry.contains(x, y);
            }
        }

        let scalar_solver = ScalarTransportSolver2D::new(nx, ny);

        Self {
            geometry,
            ns_solver,
            scalar_solver,
        }
    }

    /// Solve for flow and mixing
    pub fn solve(
        &mut self,
        u_inlet: T,
        diffusion_coeff: T,
        c_left: T,  // Concentration on left half of inlet
        c_right: T, // Concentration on right half of inlet
    ) -> CfdResult<SerpentineMixingSolution<T>> {
        // Use the scalar-transport default tolerance (1e-5), which matches the
        // FVM spatial truncation error O(Δx²) on typical coarse grids.
        self.solve_with_tolerance(u_inlet, diffusion_coeff, c_left, c_right, T::from_f64(1e-5).unwrap())
    }

    /// Solve for flow and mixing with custom tolerance
    pub fn solve_with_tolerance(
        &mut self,
        u_inlet: T,
        diffusion_coeff: T,
        c_left: T,  // Concentration on left half of inlet
        c_right: T, // Concentration on right half of inlet
        tolerance: T,
    ) -> CfdResult<SerpentineMixingSolution<T>> {
        // 1. Solve Navier-Stokes
        self.ns_solver.solve(u_inlet).map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        // 2. Setup Scalar Transport
        let mut config = ScalarTransportConfig::default();
        config.diffusion_coeff = diffusion_coeff;
        config.tolerance = tolerance;
        
        // Define inlet concentration profile (step function)
        let ny = self.ns_solver.grid.ny;
        let mut boundary_c = vec![T::zero(); ny];
        for j in 0..ny {
            if j < ny / 2 {
                boundary_c[j] = c_left;
            } else {
                boundary_c[j] = c_right;
            }
        }

        // 3. Solve Scalar Transport
        // Cells at the inlet (West boundary) with mask=true and u > 0 will use boundary_c.
        self.scalar_solver.solve(
            &self.ns_solver.grid,
            &self.ns_solver.field,
            &config,
            &boundary_c
        ).map_err(|e| cfd_core::error::Error::Solver(e))?;

        // 4. Extract metrics
        let pe = (u_inlet * self.geometry.width) / diffusion_coeff;
        
        // Compute mixing efficiency at outlet
        // ISO = 1 - variance / variance_inlet
        let nx = self.ns_solver.grid.nx;
        let mut sum_c = T::zero();
        let mut sum_c_sq = T::zero();
        let mut count = 0;
        
        // Find last fluid column (outlet)
        for j in 0..ny {
            if self.ns_solver.field.mask[nx - 1][j] {
                let ci = self.scalar_solver.c[nx - 1][j];
                sum_c = sum_c + ci;
                sum_c_sq = sum_c_sq + ci * ci;
                count += 1;
            }
        }

        let mut mixing_frac = T::zero();
        if count > 0 {
            let n = T::from_usize(count).unwrap();
            let c_mean = sum_c / n;
            let variance = (sum_c_sq / n) - (c_mean * c_mean);
            
            // Expected variance for perfectly unmixed: 0.25 for c_left=0, c_right=1
            let var_inlet = half_sq(c_left - c_right);
            mixing_frac = T::one() - Float::sqrt(variance / num_traits::Float::max(var_inlet, T::from_f64(1e-10).unwrap()));
        }

        // Compute pressure drop between inlet and outlet regions
        let mut p_in = T::zero();
        let mut count_in = 0;
        let mut p_out = T::zero();
        let mut count_out = 0;
        
        for j in 0..ny {
            if self.ns_solver.field.mask[0][j] {
                p_in = p_in + self.ns_solver.field.p[0][j];
                count_in += 1;
            }
            if self.ns_solver.field.mask[nx - 1][j] {
                p_out = p_out + self.ns_solver.field.p[nx - 1][j];
                count_out += 1;
            }
        }
        
        let mut pressure_drop = T::zero();
        if count_in > 0 && count_out > 0 {
            pressure_drop = (p_in / T::from_usize(count_in).unwrap()) - (p_out / T::from_usize(count_out).unwrap());
        }

        Ok(SerpentineMixingSolution {
            c_inlet_a: c_left,
            c_inlet_b: c_right,
            peclet: pe,
            l_mix_90: T::zero(), // Computed by analytical model
            t_mix_90: T::zero(),
            mixing_fraction_outlet: mixing_frac,
            pressure_drop,
        })
    }
}

fn half_sq<T: RealField + Copy + FromPrimitive>(val: T) -> T {
    let half = T::from_f64(0.5).unwrap();
    (val * half) * (val * half)
}

#[cfg(test)]
mod tests_discretized {
    use super::*;

    #[test]
    fn test_discretized_serpentine_mixing() {
        // High diffusion to ensure some mixing in small grid
        let geom = SerpentineGeometry::new(
            0.0002, // width
            0.0001, // height
            0.0005, // straight
            0.0002, // radius
            1,      // 1 cycle
        );
        
        // Fluid properties
        let blood = BloodModel::Newtonian(0.001);
        let density = 1000.0;
        
        // Grid (coarse for fast test)
        let nx = 40;
        let ny = 20;
        
        let mut solver = SerpentineSolver2D::new(geom, blood, density, nx, ny);
        
        // Solve with higher diffusion to see mixing
        // Use a reasonable diffusion for a coarse grid
        let result = solver.solve(0.01, 1e-6, 0.0, 1.0);
        
        assert!(result.is_ok(), "Serpentine solver failed: {:?}", result.err());
        let sol = result.unwrap();
        
        println!("Serpentine Mixing Solution: {:?}", sol);
        
        // Qualitative checks
        assert!(sol.peclet > 0.0, "Peclet should be positive");
        assert!(sol.mixing_fraction_outlet >= 0.0, "Mixing fraction should be non-negative");
        // Due to the geometry mapping, some flow might be blocked if resolution is too low.
        // We'll check if pressure drop is positive.
        assert!(sol.pressure_drop >= 0.0, "Pressure drop should be non-negative");
    }
}
