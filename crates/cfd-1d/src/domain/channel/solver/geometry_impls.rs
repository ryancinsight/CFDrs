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
use crate::scalar::Cfd1dScalar;
use aequitas::systems::si::quantities::{Area, Length};
use cfd_core::physics::constants::mathematical::{PI, numeric};
use eunomia::{FloatElement, NumericElement};

impl<T: Cfd1dScalar + Copy + FloatElement> ChannelGeometry<T> {
    /// Create a rectangular channel geometry.
    pub fn rectangular(length: T, width: T, height: T, roughness: T) -> Self {
        use crate::domain::channel::surface::{SurfaceProperties, Wettability};
        Self {
            channel_type: ChannelType::Straight,
            length: Length::from_base(length),
            cross_section: CrossSection::Rectangular {
                width: Length::from_base(width),
                height: Length::from_base(height),
            },
            surface: SurfaceProperties {
                roughness: Length::from_base(roughness),
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
            length: Length::from_base(length),
            cross_section: CrossSection::Circular {
                diameter: Length::from_base(diameter),
            },
            surface: SurfaceProperties {
                roughness: Length::from_base(roughness),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    }

    /// Cross-sectional area `A` \[m²].
    pub fn area(&self) -> Area<T> {
        let area = match &self.cross_section {
            CrossSection::Rectangular { width, height } => width.into_base() * height.into_base(),
            CrossSection::Circular { diameter } => {
                let pi = T::from_f64_or_zero(PI);
                let two = T::from_f64_or_zero(numeric::TWO);
                let radius = diameter.into_base() / two;
                pi * radius * radius
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                let pi = T::from_f64_or_zero(PI);
                let four = T::from_f64_or_zero(numeric::FOUR);
                pi * major_axis.into_base() * minor_axis.into_base() / four
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                (top_width.into_base() + bottom_width.into_base()) * height.into_base()
                    / (T::one() + T::one())
            }
            CrossSection::Custom { area, .. } => area.into_base(),
        };
        Area::from_base(area)
    }

    /// Hydraulic diameter `D_h = 4A / P` \[m].
    pub fn hydraulic_diameter(&self) -> Length<T> {
        let hydraulic_diameter = match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area().into_base()
                    / ((T::one() + T::one()) * (width.into_base() + height.into_base()))
            }
            CrossSection::Circular { diameter } => diameter.into_base(),
            CrossSection::Elliptical { .. } => {
                let four = T::one() + T::one() + T::one() + T::one();
                four * self.area().into_base() / self.wetted_perimeter().into_base()
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let area = self.area().into_base();
                let hw = (top_width.into_base() - bottom_width.into_base()) / (T::one() + T::one());
                let side_length = <T as NumericElement>::sqrt(
                    <T as FloatElement>::powi(height.into_base(), 2)
                        + <T as FloatElement>::powi(hw, 2),
                );
                let perimeter = top_width.into_base()
                    + bottom_width.into_base()
                    + (T::one() + T::one()) * side_length;
                (T::one() + T::one() + T::one() + T::one()) * area / perimeter
            }
            CrossSection::Custom {
                hydraulic_diameter, ..
            } => hydraulic_diameter.into_base(),
        };
        Length::from_base(hydraulic_diameter)
    }

    /// Wetted perimeter `P` \[m].
    ///
    /// For elliptical cross-sections, this uses the exact AGM method for the
    /// complete elliptic integral of the second kind (see module-level docs).
    pub fn wetted_perimeter(&self) -> Length<T> {
        let perimeter = match &self.cross_section {
            CrossSection::Rectangular { width, height } => {
                (T::one() + T::one()) * (width.into_base() + height.into_base())
            }
            CrossSection::Circular { diameter } => T::pi() * diameter.into_base(),
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => self.ellipse_perimeter_agm(major_axis.into_base(), minor_axis.into_base()),
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let hw = (top_width.into_base() - bottom_width.into_base()) / (T::one() + T::one());
                let side_length = <T as NumericElement>::sqrt(
                    <T as FloatElement>::powi(height.into_base(), 2)
                        + <T as FloatElement>::powi(hw, 2),
                );
                top_width.into_base()
                    + bottom_width.into_base()
                    + (T::one() + T::one()) * side_length
            }
            CrossSection::Custom {
                area,
                hydraulic_diameter,
            } => {
                (T::one() + T::one() + T::one() + T::one()) * area.into_base()
                    / hydraulic_diameter.into_base()
            }
        };
        Length::from_base(perimeter)
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
        let mut b_n = <T as NumericElement>::sqrt(T::one() - m);
        let mut c_n = <T as NumericElement>::sqrt(m);
        let mut sum = c_n * c_n / two;
        let mut power = T::one();
        let tolerance = T::from_f64_or_one(1e-14);

        for _ in 0..20 {
            let a_next = (a_n + b_n) / two;
            let b_next = <T as NumericElement>::sqrt(a_n * b_n);
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
