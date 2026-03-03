//! TPMS fill specification for shell cuboid cavities.
//!
//! [`TpmsSurfaceKind`] enumerates the nine TPMS surface types available in
//! `cfd-mesh`, and [`TpmsFillSpec`] bundles the surface choice with lattice
//! parameters (period, iso-value, resolution) needed by the downstream
//! marching-cubes extraction.
//!
//! These types are serialisable to JSON for interchange with mesh and
//! simulation pipelines.

use serde::{Deserialize, Serialize};

use crate::error::{GeometryError, GeometryResult};

// ── TPMS surface enumeration ──────────────────────────────────────────────────

/// Which TPMS surface fills the shell cavity.
///
/// Each variant maps 1-to-1 to a [`Tpms`](https://docs.rs/cfd-mesh)
/// implementor in `cfd-mesh::domain::geometry::tpms`.
///
/// | Variant | Discoverer | Implicit form |
/// |---------|------------|---------------|
/// | Gyroid | Schoen 1970 | sin(kx)cos(ky)+sin(ky)cos(kz)+sin(kz)cos(kx) |
/// | SchwarzP | Schwarz 1890 | cos(kx)+cos(ky)+cos(kz) |
/// | SchwarzD | Schwarz 1890 | sin(kx)sin(ky)sin(kz)+… |
/// | Neovius | Neovius 1883 | 3Σcos+4Πcos |
/// | Lidinoid | Lidin 1990 | ½Σsin2cos·sin−½Σcos2cos2+0.15 |
/// | Iwp | Schoen 1970 | 2Σcos·cos−Σcos2 |
/// | SplitP | Fogden & Hyde 1992 | 1.1Σsin2sincos−… |
/// | Frd | Schoen 1970 | 4Πcos−Σcos2cos2 |
/// | FischerKochCY | Fischer & Koch 1989 | Σcos2·sin·cos |
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TpmsSurfaceKind {
    /// Schoen gyroid: sin(kx)cos(ky)+sin(ky)cos(kz)+sin(kz)cos(kx).
    Gyroid,
    /// Schwarz P-surface: cos(kx)+cos(ky)+cos(kz).
    SchwarzP,
    /// Schwarz D-surface: sin(kx)sin(ky)sin(kz)+…
    SchwarzD,
    /// Neovius surface: 3(cos(kx)+cos(ky)+cos(kz))+4cos(kx)cos(ky)cos(kz).
    Neovius,
    /// Lidinoid surface.
    Lidinoid,
    /// I-WP (Schoen): 2(cos(kx)cos(ky)+…)−(cos(2kx)+…).
    Iwp,
    /// Split P (Fogden & Hyde): 1.1(sin(2kx)sin(kz)cos(ky)+…)−…
    SplitP,
    /// F-RD (Schoen): 4cos(kx)cos(ky)cos(kz)−(cos(2kx)cos(2ky)+…).
    Frd,
    /// Fischer-Koch C(Y): cos(2kx)sin(ky)cos(kz)+…
    FischerKochCY,
}

impl TpmsSurfaceKind {
    /// Short human-readable label for display in schematics.
    #[must_use]
    pub fn label(&self) -> &'static str {
        match self {
            Self::Gyroid => "Gyroid",
            Self::SchwarzP => "Schwarz P",
            Self::SchwarzD => "Schwarz D",
            Self::Neovius => "Neovius",
            Self::Lidinoid => "Lidinoid",
            Self::Iwp => "I-WP",
            Self::SplitP => "Split P",
            Self::Frd => "F-RD",
            Self::FischerKochCY => "Fischer-Koch CY",
        }
    }
}

// ── Adaptive gradient ─────────────────────────────────────────────────────────

/// Three-zone adaptive period gradient for size-based cell separation.
///
/// The TPMS period varies spatially across the shell cavity to create
/// size-dependent filtration analogous to asymmetric bifurcation/trifurcation
/// separation:
///
/// | Zone | X extent | Period | Purpose |
/// |------|----------|--------|---------|
/// | Separation | `0 → sep_frac` | `λ(y)`: periphery → center | RBC (7 µm) blocked at walls |
/// | Transition | `sep_frac → sep_frac + trans_frac` | Linear blend → uniform | Gradual pore equalisation |
/// | Remerging | remainder | Uniform `period_center_mm` | Blood reconverges before outlet |
///
/// ## Cell-size anchoring
///
/// - `period_periphery_mm` → pore throat < 7 µm (blocks RBCs + lymphocytes)
/// - `period_center_mm` → pore throat > 17.5 µm (passes WBCs + CTCs)
///
/// ## Invariant
///
/// `separation_fraction + transition_fraction ≤ 1.0`; the remainder is the
/// remerging zone.  Both fractions must be in `(0, 1)`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdaptiveGradient {
    /// Period at the peripheral walls [mm] (fine pores → block RBCs).
    ///
    /// Must be positive and finite.  Typical: 0.5–2.0 mm for millifluidic
    /// devices targeting 7 µm RBC exclusion.
    pub period_periphery_mm: f64,
    /// Period at the center [mm] (coarse pores → pass WBCs/CTCs).
    ///
    /// Must be positive, finite, and ≥ `period_periphery_mm`.
    /// Typical: 5.0–15.0 mm.
    pub period_center_mm: f64,
    /// Fraction of cavity X-length devoted to the graded separation zone.
    ///
    /// Range: `(0, 1)`.  Typical: 0.5–0.7.
    pub separation_fraction: f64,
    /// Fraction of cavity X-length for the linear transition zone.
    ///
    /// Range: `(0, 1)`.  `separation_fraction + transition_fraction ≤ 1.0`.
    /// Typical: 0.1–0.2.
    pub transition_fraction: f64,
}

impl AdaptiveGradient {
    /// Fraction of cavity X-length for the uniform remerging zone.
    ///
    /// `1.0 - separation_fraction - transition_fraction`.
    #[must_use]
    pub fn remerging_fraction(&self) -> f64 {
        1.0 - self.separation_fraction - self.transition_fraction
    }

    /// Compute the local TPMS period at a normalised position `(x_frac, y_frac)`.
    ///
    /// - `x_frac ∈ [0, 1]` — fractional position along the cavity X-axis
    ///   (0 = inlet side, 1 = outlet side).
    /// - `y_frac ∈ [0, 1]` — fractional position across the cavity width
    ///   (0 = bottom wall, 1 = top wall, 0.5 = center).
    ///
    /// Returns the period [mm] at this spatial position.
    #[must_use]
    pub fn period_at(&self, x_frac: f64, y_frac: f64) -> f64 {
        let p_periph = self.period_periphery_mm;
        let p_center = self.period_center_mm;
        let sf = self.separation_fraction;
        let tf = self.transition_fraction;

        // Y-gradient: distance from center (0 at center, 1 at walls).
        let y_wall_dist = (2.0 * (y_frac - 0.5)).abs(); // 0 at center, 1 at walls
        let y_period = p_center + (p_periph - p_center) * y_wall_dist;

        if x_frac <= sf {
            // Separation zone: full Y-gradient.
            y_period
        } else if x_frac <= sf + tf {
            // Transition zone: linear blend from Y-gradient → uniform center.
            let blend = (x_frac - sf) / tf;
            y_period * (1.0 - blend) + p_center * blend
        } else {
            // Remerging zone: uniform center period.
            p_center
        }
    }

    /// Validate all gradient parameters.
    pub(crate) fn validate(&self) -> GeometryResult<()> {
        if !self.period_periphery_mm.is_finite() || self.period_periphery_mm <= 0.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "period_periphery_mm must be finite and positive, got {}",
                    self.period_periphery_mm
                ),
            });
        }
        if !self.period_center_mm.is_finite() || self.period_center_mm <= 0.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "period_center_mm must be finite and positive, got {}",
                    self.period_center_mm
                ),
            });
        }
        if self.period_center_mm < self.period_periphery_mm {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "period_center_mm ({}) must be >= period_periphery_mm ({})",
                    self.period_center_mm, self.period_periphery_mm
                ),
            });
        }
        let sf = self.separation_fraction;
        let tf = self.transition_fraction;
        if !sf.is_finite() || sf <= 0.0 || sf >= 1.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!("separation_fraction must be in (0, 1), got {sf}"),
            });
        }
        if !tf.is_finite() || tf <= 0.0 || tf >= 1.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!("transition_fraction must be in (0, 1), got {tf}"),
            });
        }
        if sf + tf > 1.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "separation_fraction ({sf}) + transition_fraction ({tf}) = {} > 1.0",
                    sf + tf
                ),
            });
        }
        Ok(())
    }
}

impl Default for AdaptiveGradient {
    /// Default: 60% separation, 20% transition, 20% remerging.
    /// Periphery 1.5 mm (fine pores), center 8.0 mm (coarse pores).
    fn default() -> Self {
        Self {
            period_periphery_mm: 1.5,
            period_center_mm: 8.0,
            separation_fraction: 0.6,
            transition_fraction: 0.2,
        }
    }
}

// ── Fill specification ────────────────────────────────────────────────────────

/// Parameters specifying how a TPMS surface fills a shell cavity.
///
/// This struct is designed for JSON interchange: `cfd-schematics` produces it,
/// and `cfd-mesh` consumes it via the `ShellMeshPipeline`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TpmsFillSpec {
    /// Which TPMS surface to extract.
    pub surface: TpmsSurfaceKind,
    /// Baseline unit-cell period [mm] (used when `gradient` is `None`).
    ///
    /// Typical range for millifluidic devices: 2–10 mm.
    pub period_mm: f64,
    /// Level-set iso-value.  `0.0` = the exact minimal surface mid-sheet.
    /// Positive values thicken the lattice walls; negative values thin them.
    pub iso_value: f64,
    /// Marching-cubes voxels per axis.  Higher = denser mesh, longer build.
    /// Must be ≥ 4.  Typical: 32–128.
    pub resolution: usize,
    /// Optional adaptive period gradient for size-based cell separation.
    ///
    /// When `Some`, the period varies spatially across the cavity
    /// (fine at periphery, coarse at center, uniform at outlet).
    /// When `None`, the uniform `period_mm` is used everywhere.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gradient: Option<AdaptiveGradient>,
}

impl Default for TpmsFillSpec {
    /// Gyroid at 5 mm period, mid-surface, resolution 64.
    fn default() -> Self {
        Self {
            surface: TpmsSurfaceKind::Gyroid,
            period_mm: 5.0,
            iso_value: 0.0,
            resolution: 64,
            gradient: None,
        }
    }
}

impl TpmsFillSpec {
    /// Validate parameter ranges.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError`] if period is non-positive / non-finite
    /// or resolution is below 4.
    pub fn validate(&self) -> GeometryResult<()> {
        if !self.period_mm.is_finite() || self.period_mm <= 0.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "tpms_fill period_mm must be a finite positive value, got {}",
                    self.period_mm
                ),
            });
        }
        if self.resolution < 4 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!("tpms_fill resolution must be >= 4, got {}", self.resolution),
            });
        }
        if !self.iso_value.is_finite() {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!("tpms_fill iso_value must be finite, got {}", self.iso_value),
            });
        }
        if let Some(ref grad) = self.gradient {
            grad.validate()?;
        }
        Ok(())
    }

    /// Whether this fill uses an adaptive gradient vs uniform period.
    #[must_use]
    pub fn is_graded(&self) -> bool {
        self.gradient.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_fill_spec_validates() {
        TpmsFillSpec::default()
            .validate()
            .expect("default fill spec should be valid");
    }

    #[test]
    fn fill_spec_json_roundtrip() {
        let spec = TpmsFillSpec {
            surface: TpmsSurfaceKind::SchwarzP,
            period_mm: 3.0,
            iso_value: 0.1,
            resolution: 32,
            gradient: None,
        };
        let json = serde_json::to_string(&spec).unwrap();
        let parsed: TpmsFillSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, parsed);
    }

    #[test]
    fn fill_spec_rejects_zero_period() {
        let spec = TpmsFillSpec {
            period_mm: 0.0,
            ..TpmsFillSpec::default()
        };
        assert!(spec.validate().is_err());
    }

    #[test]
    fn fill_spec_rejects_low_resolution() {
        let spec = TpmsFillSpec {
            resolution: 2,
            ..TpmsFillSpec::default()
        };
        assert!(spec.validate().is_err());
    }

    #[test]
    fn surface_kind_labels_nonempty() {
        let kinds = [
            TpmsSurfaceKind::Gyroid,
            TpmsSurfaceKind::SchwarzP,
            TpmsSurfaceKind::SchwarzD,
            TpmsSurfaceKind::Neovius,
            TpmsSurfaceKind::Lidinoid,
            TpmsSurfaceKind::Iwp,
            TpmsSurfaceKind::SplitP,
            TpmsSurfaceKind::Frd,
            TpmsSurfaceKind::FischerKochCY,
        ];
        for k in &kinds {
            assert!(!k.label().is_empty());
        }
    }

    // ── Adaptive gradient tests ───────────────────────────────────────────

    fn make_default_gradient() -> AdaptiveGradient {
        AdaptiveGradient::default()
    }

    #[test]
    fn default_gradient_validates() {
        make_default_gradient()
            .validate()
            .expect("default gradient should be valid");
    }

    #[test]
    fn gradient_remerging_fraction() {
        let g = make_default_gradient();
        let rem = g.remerging_fraction();
        assert!(
            (rem - 0.2).abs() < 1e-9,
            "default: 1.0 - 0.6 - 0.2 = 0.2, got {rem}"
        );
    }

    #[test]
    fn gradient_period_at_separation_center() {
        let g = make_default_gradient();
        // (x=0.3 → separation zone, y=0.5 → center) → period_center_mm
        let p = g.period_at(0.3, 0.5);
        assert!(
            (p - g.period_center_mm).abs() < 1e-9,
            "center of separation zone should be center period"
        );
    }

    #[test]
    fn gradient_period_at_separation_wall() {
        let g = make_default_gradient();
        // (x=0.3 → separation zone, y=0.0 → bottom wall) → period_periphery_mm
        let p = g.period_at(0.3, 0.0);
        assert!(
            (p - g.period_periphery_mm).abs() < 1e-9,
            "wall of separation zone should be periphery period"
        );
    }

    #[test]
    fn gradient_period_at_remerging() {
        let g = make_default_gradient();
        // x=0.9 → remerging zone → always period_center_mm regardless of y
        let p_center = g.period_at(0.9, 0.5);
        let p_wall = g.period_at(0.9, 0.0);
        assert!((p_center - g.period_center_mm).abs() < 1e-9);
        assert!((p_wall - g.period_center_mm).abs() < 1e-9);
    }

    #[test]
    fn gradient_transition_blends() {
        let g = make_default_gradient();
        // x=0.7 → midway through transition (0.6..0.8), y=0.0 → wall
        let p = g.period_at(0.7, 0.0);
        // Expected: 50% blend between periphery (at wall) and center
        let expected = g.period_periphery_mm * 0.5 + g.period_center_mm * 0.5;
        assert!(
            (p - expected).abs() < 1e-9,
            "midpoint transition should be 50/50 blend, got {p}"
        );
    }

    #[test]
    fn gradient_rejects_zone_sum_over_one() {
        let g = AdaptiveGradient {
            separation_fraction: 0.7,
            transition_fraction: 0.4,
            ..AdaptiveGradient::default()
        };
        assert!(g.validate().is_err(), "sep + trans > 1.0 must fail");
    }

    #[test]
    fn gradient_rejects_center_less_than_periphery() {
        let g = AdaptiveGradient {
            period_periphery_mm: 10.0,
            period_center_mm: 5.0,
            ..AdaptiveGradient::default()
        };
        assert!(g.validate().is_err(), "center < periphery must fail");
    }

    #[test]
    fn fill_spec_with_gradient_json_roundtrip() {
        let spec = TpmsFillSpec {
            gradient: Some(AdaptiveGradient::default()),
            ..TpmsFillSpec::default()
        };
        let json = serde_json::to_string_pretty(&spec).unwrap();
        let parsed: TpmsFillSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, parsed);
        assert!(parsed.is_graded());
    }

    #[test]
    fn fill_spec_without_gradient_not_graded() {
        assert!(!TpmsFillSpec::default().is_graded());
    }
}
