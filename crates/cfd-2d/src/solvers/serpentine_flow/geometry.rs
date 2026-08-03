use crate::scalar::Cfd2dScalar;
use crate::scalar::{self, from_f64};
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Serpentine channel geometry with periodic turns
///
/// Defines a channel that alternates between straight sections and 90° turns.
/// The path creates a snake-like pattern that enhances mixing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineGeometry<T: Cfd2dScalar + Copy> {
    /// Channel width \[m]
    pub width: T,
    /// Channel height \[m] (constant)
    pub height: T,
    /// Straight section length \[m]
    pub l_straight: T,
    /// Turn radius \[m] (for curved turns)
    pub turn_radius: T,
    /// Number of complete back-and-forth cycles
    pub n_cycles: usize,
}

impl<T: Cfd2dScalar + Copy + FloatElement> SerpentineGeometry<T> {
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
            width: from_f64::<T>(200e-6),
            height: from_f64::<T>(50e-6),
            l_straight: from_f64::<T>(500e-6),
            turn_radius: from_f64::<T>(200e-6),
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
        let pi = from_f64::<T>(PI);
        let two = from_f64::<T>(2.0);
        let turn_length = pi / from_f64::<T>(2.0) * self.turn_radius;

        scalar::from_usize::<T>(self.n_cycles) * (two * self.l_straight + turn_length)
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
        let three = from_f64::<T>(3.0);
        (self.l_straight / (self.width + from_f64::<T>(1e-15))) * three * peclet
    }

    /// Get bounding box [min_x, max_x, min_y, max_y]
    pub fn bounding_box(&self) -> [T; 4] {
        let r = self.turn_radius;
        let ls = self.l_straight;
        let w = self.width;
        let half_w = w / from_f64::<T>(2.0);
        let four_r = r * from_f64::<T>(4.0);

        [
            -r - half_w,
            ls + r + half_w,
            scalar::zero(),
            four_r * scalar::from_usize::<T>(self.n_cycles),
        ]
    }
}

impl<T: Cfd2dScalar + Copy + FloatElement + std::ops::Rem<Output = T>> SerpentineGeometry<T> {
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
        let half_w = w / from_f64::<T>(2.0);
        let four_r = r * from_f64::<T>(4.0);

        let y_in_cycle = y % four_r;

        let inner_r = r - half_w;
        let outer_r = r + half_w;
        let inner_r_sq = inner_r * inner_r;
        let outer_r_sq = outer_r * outer_r;

        if y_in_cycle < r + half_w && y_in_cycle > r - half_w && x >= -r - half_w && x <= ls {
            return true;
        }

        if y_in_cycle < from_f64::<T>(3.0) * r + half_w
            && y_in_cycle > from_f64::<T>(3.0) * r - half_w
            && x >= scalar::zero()
            && x <= ls
        {
            return true;
        }

        let dx_r = x - ls;
        let dy_r = y_in_cycle - r * from_f64::<T>(2.0);
        let d_sq_r = dx_r * dx_r + dy_r * dy_r;
        if x >= ls && d_sq_r >= inner_r_sq && d_sq_r <= outer_r_sq {
            return true;
        }

        let dy_l_top = y_in_cycle - four_r;
        let d_sq_l_top = x * x + dy_l_top * dy_l_top;
        if x <= scalar::zero() && d_sq_l_top >= inner_r_sq && d_sq_l_top <= outer_r_sq {
            return true;
        }

        let d_sq_l_bot = x * x + y_in_cycle * y_in_cycle;
        if x <= scalar::zero() && d_sq_l_bot >= inner_r_sq && d_sq_l_bot <= outer_r_sq {
            return true;
        }

        false
    }
}
