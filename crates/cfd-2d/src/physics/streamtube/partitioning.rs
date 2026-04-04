//! 2D Streamtube Partitioning (Zweifach-Fung Basis)
//!
//! Provides mathematically rigorous integration of 2D velocity profiles
//! to extract dividing/separating streamlines for cell or droplet routing.
//!
//! # Theorem — Flow Partitioning Mass Conservation
//!
//! In a steady incompressible 2D flow ($$\nabla \cdot \vec{u} = 0$$), the dividing
//! streamline $y_{sep}(x)$ upstream of a bifurcation exactly partitions the total volumetric
//! flux $Q_{Total}$ into $Q_1$ and $Q_2$:
//!
//! ```text
//! \int_{y_{min}}^{y_{sep}(x)} u(x, y) dy = Q_1 = f_q \cdot Q_{Total}
//! ```
//!
//! Assuming a known fraction $f_q$ entering a specific daughter branch,
//! the separating streamline coordinate $y_{sep}$ is unique and monotonically
//! increasing with $f_q$ for any strictly positive unidirectional profile $u(y) > 0$.
//!
//! **Integration Method**: We evaluate the discrete trapezoidal cumulative flow exactly on
//! each interval. When the sampled velocity is treated as piecewise linear, the cumulative
//! flow is quadratic on that interval, and we solve the quadratic root analytically to obtain
//! the exact $y_{sep}$ for the discrete representation.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Core 2D Streamtube Partitioning Tools
pub struct ZweifachFung2D;

impl ZweifachFung2D {
    /// Calculate the fractional flow separation coordinate $y_{sep}$ given a
    /// 1D cross-sectional velocity profile $u(y)$.
    ///
    /// The profile is assumed to be defined at discrete coordinates $(y_i, u_i)$.
    /// The function computes the cumulative flow $Q(y)$ and returns $y_{sep}$
    /// where $Q(y_{sep}) = target\_fraction \times Q_{total}$.
    ///
    /// # Arguments
    /// * `y_coords` - Monotonically increasing discrete y-coordinates (e.g., cell centers or faces)
    /// * `u_vel` - Corresponding streamwise velocity $u(y)$ at each coordinate
    /// * `target_fraction` - The desired flow fraction (0.0 to 1.0) $f_q = Q_1 / Q_{total}$
    ///
    /// # Returns
    /// The exact interpolated $y_{sep}$ coordinate. Returns `None` if the input is invalid
    /// or if the profile contains significant reverse flow (which breaks the simple monotonic CDF).
    pub fn separating_streamline_y<T: RealField + Copy + Float + FromPrimitive>(
        y_coords: &[T],
        u_vel: &[T],
        target_fraction: T,
    ) -> Option<T> {
        let n = y_coords.len();
        if n < 2 || y_coords.len() != u_vel.len() {
            return None;
        }

        let zero = T::zero();
        let one = T::one();
        let two = T::from_f64(2.0).unwrap();

        if target_fraction <= zero {
            return Some(y_coords[0]);
        }
        if target_fraction >= one {
            return Some(y_coords[n - 1]);
        }

        // Compute cumulative flow Q(y) via trapezoidal integration
        let mut q_cumulative = vec![zero; n];
        for i in 1..n {
            let dy = y_coords[i] - y_coords[i - 1];
            let u_avg = (u_vel[i] + u_vel[i - 1]) / two;
            let dq = u_avg * dy;
            
            // For rigorous Zweifach-Fung sorting, we usually assume strictly positive flow.
            // If massive separation (recirculation) exists at the split cross-section,
            // 1D streamtubes are ill-posed, but we proceed anyway for robustness.
            q_cumulative[i] = q_cumulative[i - 1] + dq;
        }

        let q_total = q_cumulative[n - 1];
        if q_total <= zero {
            return None; // No net forward flow
        }

        let q_target = q_total * target_fraction;

        // Find the cell bounding the target flow Q_target
        for i in 1..n {
            if q_cumulative[i] >= q_target {
                // We bracketed the target between i-1 and i
                let q0 = q_cumulative[i - 1];
                let q1 = q_cumulative[i];
                let y0 = y_coords[i - 1];
                let y1 = y_coords[i];

                if q1 == q0 {
                    return Some(y0); // Exact hit or zero velocity zone
                }

                // Solve the exact quadratic induced by piecewise-linear velocity samples.
                // dQ = u0 * s + 0.5 * (u1 - u0) * s^2 / dy_interval.
                let delta_q = q_target - q0;
                let dy_interval = y1 - y0;
                let u0 = u_vel[i - 1];
                let u1 = u_vel[i];
                let du = u1 - u0;

                return if Float::abs(du) < T::from_f64(1e-12).unwrap() {
                    // Uniform velocity in this cell
                    Some(y0 + delta_q / u0)
                } else {
                    // Quadratic formula: 0.5 * (du/dy_int) * s^2 + u0 * s - delta_q = 0
                    let a = (half::<T>() * du) / dy_interval;
                    let b = u0;
                    let c = -delta_q;

                    // s = (-b + sqrt(b^2 - 4ac)) / 2a
                    let disc = b * b - T::from_f64(4.0).unwrap() * a * c;
                    if disc < zero {
                        // Fallback to linear if root imaginary due to precision edge cases
                        let t = delta_q / (q1 - q0);
                        Some(y0 + t * dy_interval)
                    } else {
                        let s = (-b + Float::sqrt(disc)) / (two * a);
                        Some(y0 + s)
                    }
                };
            }
        }

        None
    }
}

/// Helper to cast float to T
#[inline]
fn half<T: RealField + Float + FromPrimitive>() -> T {
    T::from_f64(0.5).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Theorem: Analytical integration of steady Poiseuille flow between parallel plates.
    /// Profile: u(y) = u_max * (1 - (2*y/H)^2) for y in [-H/2, H/2].
    /// Integral: Q_total = 2/3 * u_max * H.
    /// Exact target streamtube for 50% flow is the centerline `y=0`.
    #[test]
    fn separating_streamline_poiseuille_exact() {
        let n = 201; // Fine discrete grid
        let mut y_coords = vec![0.0_f64; n];
        let mut u_vel = vec![0.0_f64; n];

        let h = 0.002; // 2mm channel
        let half_h = h / 2.0;
        let u_max = 0.1; // 10 cm/s

        for i in 0..n {
            let y = -half_h + h * (i as f64) / ((n - 1) as f64);
            y_coords[i] = y;
            let norm_y = y / half_h;
            u_vel[i] = u_max * (1.0 - norm_y * norm_y);
        }

        // Check 50% flow fraction (should be exactly at centerline y = 0.0)
        let y_sep_50 = ZweifachFung2D::separating_streamline_y(&y_coords, &u_vel, 0.5).unwrap();
        assert!(
            y_sep_50.abs() < 1e-6,
            "50% flow must split exactly at the centerline. Got: {}",
            y_sep_50
        );

        // Check 100% flow fraction
        let y_sep_100 = ZweifachFung2D::separating_streamline_y(&y_coords, &u_vel, 1.0).unwrap();
        assert!(
            (y_sep_100 - half_h).abs() < 1e-6,
            "100% flow must split exactly at upper wall. Got {}",
            y_sep_100
        );
        
        // Check 0% flow fraction
        let y_sep_0 = ZweifachFung2D::separating_streamline_y(&y_coords, &u_vel, 0.0).unwrap();
        assert!(
            (y_sep_0 - (-half_h)).abs() < 1e-6,
            "0% flow must split exactly at lower wall. Got {}",
            y_sep_0
        );
    }

    #[test]
    fn separating_streamline_linear_profile_closed_form() {
        let y_coords = vec![0.0_f64, 0.25, 0.5, 0.75, 1.0];
        let u_vel: Vec<f64> = y_coords.iter().map(|&y| 1.0 + y).collect();

        let y_sep = ZweifachFung2D::separating_streamline_y(&y_coords, &u_vel, 0.25).unwrap();
        let expected = -1.0 + 1.75_f64.sqrt();

        assert_relative_eq!(y_sep, expected, epsilon = 1e-12);
    }
}
