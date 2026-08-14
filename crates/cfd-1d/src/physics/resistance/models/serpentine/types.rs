use cfd_core::CfdScalar;
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Bend type in the serpentine channel
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BendType {
    /// Sharp 180° U-turn (no fillet radius)
    Sharp,
    /// Smooth 180° bend with specified R/D ratio
    Smooth {
        /// Ratio of bend radius to hydraulic diameter (R/D_h)
        radius_to_dh_ratio: f64,
    },
}

impl BendType {
    /// Minor loss coefficient for a single bend: K = C1 + C2/Re
    ///
    /// Returns (C1, C2) constants based on Idelchik (2007) §6.2
    #[must_use]
    pub fn loss_constants(&self) -> (f64, f64) {
        match self {
            Self::Sharp => (2.2, 250.0),
            Self::Smooth { radius_to_dh_ratio } => {
                let r = *radius_to_dh_ratio;
                if r <= 1.5 {
                    (0.9, 200.0)
                } else if r <= 3.0 {
                    (0.4, 100.0)
                } else if r <= 5.0 {
                    (0.3, 75.0)
                } else {
                    (0.2, 50.0)
                }
            }
        }
    }

    /// Compute bend loss coefficient K at the given Reynolds number
    pub fn loss_coefficient<T: CfdScalar>(&self, reynolds: T) -> T {
        let (c1, c2) = self.loss_constants();
        let c1_t = <T as FloatElement>::from_f64(c1);
        let c2_t = <T as FloatElement>::from_f64(c2);
        c1_t + c2_t / reynolds
    }
}

/// Cross-section type for the serpentine channel
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SerpentineCrossSection {
    /// Circular cross-section with given diameter \[m]
    Circular {
        /// Channel diameter \[m]
        diameter: f64,
    },
    /// Rectangular cross-section with width × height
    Rectangular {
        /// Channel width \[m]
        width: f64,
        /// Channel height (depth) \[m]
        height: f64,
    },
}

impl SerpentineCrossSection {
    /// Hydraulic diameter D_h
    #[must_use]
    pub fn hydraulic_diameter(&self) -> f64 {
        match self {
            Self::Circular { diameter } => *diameter,
            Self::Rectangular { width, height } => 2.0 * width * height / (width + height),
        }
    }

    /// Cross-sectional area \[m²]
    #[must_use]
    pub fn area(&self) -> f64 {
        match self {
            Self::Circular { diameter } => std::f64::consts::PI * diameter * diameter / 4.0,
            Self::Rectangular { width, height } => width * height,
        }
    }

    /// Aspect ratio for rectangular channels (max/min dimension)
    #[must_use]
    pub fn aspect_ratio(&self) -> f64 {
        match self {
            Self::Circular { .. } => 1.0,
            Self::Rectangular { width, height } => {
                let (a, b) = if width > height {
                    (*width, *height)
                } else {
                    (*height, *width)
                };
                a / b
            }
        }
    }

    /// Shah-London correction factor for rectangular channels
    ///
    /// f·Re product for fully developed laminar flow in rectangular ducts
    /// compared to circular pipes:
    ///
    /// f·Re = C(α) where α is the aspect ratio
    ///
    /// Reference: Shah & London (1978), Table 43
    #[must_use]
    pub fn shah_london_fre_factor(&self) -> f64 {
        match self {
            Self::Circular { .. } => 1.0,
            Self::Rectangular { width, height } => {
                let alpha = {
                    let (a, b) = if width > height {
                        (*width, *height)
                    } else {
                        (*height, *width)
                    };
                    b / a
                };
                let fre_rect = 96.0
                    * (1.0 - 1.3553 * alpha + 1.9467 * alpha.powi(2) - 1.7012 * alpha.powi(3)
                        + 0.9564 * alpha.powi(4)
                        - 0.2537 * alpha.powi(5));
                fre_rect / 64.0
            }
        }
    }
}
