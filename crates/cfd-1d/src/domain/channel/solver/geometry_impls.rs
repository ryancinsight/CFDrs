//! `ChannelGeometry` cross-section computation methods.
//!
//! # Overview
//!
//! Implements the three fundamental geometric quantities needed for hydraulic
//! resistance computation on `ChannelGeometry<T>`:
//!
//! | Method              | Definition                     |
//! |---------------------|-------------------------------|
//! | `area()`            | Cross-sectional area `A`       |
//! | `hydraulic_diameter()` | `D_h = 4A / P`             |
//! | `wetted_perimeter()` | Wetted perimeter `P`          |
//!
//! Each method dispatches on the [`CrossSection`] variant, supporting circular,
//! rectangular, elliptical, trapezoidal, and custom cross-sections.
//!
//! # Theorem: Exact Ellipse Perimeter via AGM
//!
//! The perimeter of an ellipse with semi-axes `a ≥ b` is computed exactly
//! (to machine precision) using the Arithmetic-Geometric Mean (AGM) method
//! for the complete elliptic integral of the second kind `E(m)`:
//!
//! ```text
//! P = 4a · E(m),   m = 1 − (b/a)²
//! ```
//!
//! The AGM iteration converges quadratically, reaching `f64` precision in
//! ≤ 20 iterations for all `m ∈ [0, 1)`.
//!
//! **Proof sketch**: The AGM iteration doubles the number of correct digits
//! at each step. Since `|c_0| ≤ 1`, after `n` iterations
//! `|c_n| ≤ 2^{-2^n}`, which is below `f64` epsilon after ~5 iterations.
//! We use 20 iterations for a conservative bound.
//!
//! # References
//! - Borwein, J. M. & Borwein, P. B. (1987). *Pi and the AGM*.

use crate::domain::channel::cross_section::CrossSection;
use crate::domain::channel::geometry::{ChannelGeometry, ChannelType};
use cfd_core::conversion::SafeFromF64;
use cfd_core::physics::constants::mathematical::{numeric, PI};
use nalgebra::RealField;
use num_traits::{cast::FromPrimitive, Float};

impl<T: RealField + Copy + FromPrimitive + Float> ChannelGeometry<T> {
    /// Create a rectangular channel geometry.
    pub fn rectangular(length: T, width: T, height: T, roughness: T) -> Self {
        use crate::domain::channel::surface::{SurfaceProperties, Wettability};
        Self {
            channel_type: ChannelType::Straight,
            length,
            cross_section: CrossSection::Rectangular { width, height },
            surface: SurfaceProperties {
                roughness,
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    }

    /// Create a circular channel geometry.
    pub fn circular(length: T, diameter: T, roughness: T) -> Self {
        use crate::domain::channel::surface::{SurfaceProperties, Wettability};
        Self {
            channel_type: ChannelType::Straight,
            length,
            cross_section: CrossSection::Circular { diameter },
            surface: SurfaceProperties {
                roughness,
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    }

    /// Cross-sectional area `A` [m²].
    pub fn area(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => *width * *height,
            CrossSection::Circular { diameter } => {
                let pi = T::from_f64_or_zero(PI);
                let two = T::from_f64_or_zero(numeric::TWO);
                let radius = *diameter / two;
                pi * radius * radius
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                let pi = T::from_f64_or_zero(PI);
                let four = T::from_f64_or_zero(numeric::FOUR);
                pi * *major_axis * *minor_axis / four
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => (*top_width + *bottom_width) * *height / (T::one() + T::one()),
            CrossSection::Custom { area, .. } => *area,
        }
    }

    /// Hydraulic diameter `D_h = 4A / P` [m].
    pub fn hydraulic_diameter(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area() / ((T::one() + T::one()) * (*width + *height))
            }
            CrossSection::Circular { diameter } => *diameter,
            CrossSection::Elliptical { .. } => {
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area() / self.wetted_perimeter()
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let area = self.area();
                let hw = (*top_width - *bottom_width) / (T::one() + T::one());
                let side_length = Float::sqrt(Float::powi(*height, 2) + Float::powi(hw, 2));
                let perimeter = *top_width + *bottom_width + (T::one() + T::one()) * side_length;
                (T::one() + T::one() + T::one() + T::one()) * area / perimeter
            }
            CrossSection::Custom {
                hydraulic_diameter, ..
            } => *hydraulic_diameter,
        }
    }

    /// Wetted perimeter `P` [m].
    ///
    /// For elliptical cross-sections, this uses the exact AGM method for the
    /// complete elliptic integral of the second kind (see module-level docs).
    pub fn wetted_perimeter(&self) -> T {
        match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                (T::one() + T::one()) * (*width + *height)
            }
            CrossSection::Circular { diameter } => T::pi() * *diameter,
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => self.ellipse_perimeter_agm(*major_axis, *minor_axis),
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let hw = (*top_width - *bottom_width) / (T::one() + T::one());
                let side_length = Float::sqrt(Float::powi(*height, 2) + Float::powi(hw, 2));
                *top_width + *bottom_width + (T::one() + T::one()) * side_length
            }
            CrossSection::Custom {
                area,
                hydraulic_diameter,
            } => (T::one() + T::one() + T::one() + T::one()) * *area / *hydraulic_diameter,
        }
    }

    /// Exact ellipse perimeter via the Arithmetic-Geometric Mean (AGM) method
    /// for the complete elliptic integral of the second kind.
    ///
    /// `P = 4a · E(m)`, where `m = 1 − (b/a)²`.
    fn ellipse_perimeter_agm(&self, major_axis: T, minor_axis: T) -> T {
        let pi = T::pi();
        let two = T::one() + T::one();
        let a_val = major_axis / two;
        let b_val = minor_axis / two;

        let (a, b) = if a_val > b_val {
            (a_val, b_val)
        } else {
            (b_val, a_val)
        };

        if a == b || b == T::zero() {
            return two * pi * a;
        }

        let m = T::one() - (b * b) / (a * a);

        let mut a_n = T::one();
        let mut b_n = Float::sqrt(T::one() - m);
        let mut c_n = Float::sqrt(m);
        let mut sum = c_n * c_n / two;
        let mut power = T::one();
        let tolerance = T::from_f64(1e-14).expect("AGM tolerance constant");

        for _ in 0..20 {
            let a_next = (a_n + b_n) / two;
            let b_next = Float::sqrt(a_n * b_n);
            let c_next = (a_n - b_n) / two;

            a_n = a_next;
            b_n = b_next;
            c_n = c_next;

            sum += power * c_n * c_n;
            power *= two;

            if c_n < tolerance || c_n == T::zero() {
                break;
            }
        }

        let e_m = (pi / (two * a_n)) * (T::one() - sum);
        let four = T::one() + T::one() + T::one() + T::one();
        four * a * e_m
    }
}
