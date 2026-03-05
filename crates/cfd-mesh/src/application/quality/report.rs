//! Extended mesh quality report with per-metric histograms and aggregate counts.
//!
//! [`FullQualityReport`] extends the existing [`QualityReport`] with:
//! - Fixed-width histograms for edge length, minimum angle, and aspect ratio.
//! - `bad_angle_count` — faces with min-angle below the configured threshold.
//! - `bad_aspect_count` — faces with aspect ratio above the configured threshold.
//!
//! The base [`QualityReport`] (from `validation.rs`) remains unchanged and is
//! embedded here verbatim so callers that already use `QualityReport` can
//! upgrade by accessing the `.base` field.

use crate::application::quality::histograms::Histogram;
use crate::application::quality::validation::QualityReport;
use crate::domain::core::scalar::Real;

// ── Full quality report ────────────────────────────────────────────────────────

/// Extended mesh quality report.
///
/// Contains the base per-metric statistics, per-metric histograms, and
/// aggregate counts of faces that violate standard CFD quality thresholds.
#[derive(Clone, Debug)]
pub struct FullQualityReport {
    /// Base statistics (min/max/mean for each metric).
    pub base: QualityReport,

    /// Edge-length distribution histogram.
    ///
    /// Aggregates all three edge lengths for every face.
    pub edge_length_histogram: Option<Histogram>,

    /// Minimum interior angle distribution (radians).
    pub min_angle_histogram: Option<Histogram>,

    /// Aspect-ratio distribution.
    pub aspect_ratio_histogram: Option<Histogram>,

    /// Faces whose minimum angle is below the threshold.
    pub bad_angle_count: usize,

    /// Faces whose aspect ratio exceeds the threshold.
    pub bad_aspect_count: usize,
}

impl FullQualityReport {
    /// Fraction of faces failing the minimum-angle threshold `[0, 1]`.
    #[must_use]
    pub fn bad_angle_fraction(&self) -> Real {
        if self.base.total_faces == 0 {
            return 0.0;
        }
        self.bad_angle_count as Real / self.base.total_faces as Real
    }

    /// Fraction of faces failing the aspect-ratio threshold `[0, 1]`.
    #[must_use]
    pub fn bad_aspect_fraction(&self) -> Real {
        if self.base.total_faces == 0 {
            return 0.0;
        }
        self.bad_aspect_count as Real / self.base.total_faces as Real
    }

    /// `true` when all faces satisfy both the minimum-angle and aspect-ratio thresholds.
    #[must_use]
    pub fn is_fully_acceptable(&self) -> bool {
        self.bad_angle_count == 0 && self.bad_aspect_count == 0
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::quality::validation::QualityReport;

    fn make_report(bad_angle: usize, bad_aspect: usize, total: usize) -> FullQualityReport {
        FullQualityReport {
            base: QualityReport {
                aspect_ratio: None,
                min_angle: None,
                skewness: None,
                edge_ratio: None,
                failing_faces: bad_angle + bad_aspect,
                total_faces: total,
                passed: bad_angle + bad_aspect == 0,
            },
            edge_length_histogram: None,
            min_angle_histogram: None,
            aspect_ratio_histogram: None,
            bad_angle_count: bad_angle,
            bad_aspect_count: bad_aspect,
        }
    }

    #[test]
    fn bad_fractions_correct() {
        let r = make_report(10, 5, 100);
        assert!((r.bad_angle_fraction() - 0.10).abs() < 1e-10);
        assert!((r.bad_aspect_fraction() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn is_fully_acceptable_when_no_bad_faces() {
        let r = make_report(0, 0, 50);
        assert!(r.is_fully_acceptable());
    }

    #[test]
    fn is_not_acceptable_with_bad_faces() {
        let r = make_report(1, 0, 50);
        assert!(!r.is_fully_acceptable());
    }

    #[test]
    fn bad_fraction_zero_for_empty_mesh() {
        let r = make_report(0, 0, 0);
        assert_eq!(r.bad_angle_fraction(), 0.0);
        assert_eq!(r.bad_aspect_fraction(), 0.0);
    }
}
