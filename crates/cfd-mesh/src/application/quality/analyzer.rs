//! DIP-compliant quality analysis trait and standard implementation.
//!
//! ## Design
//!
//! [`QualityAnalyzer`] is a trait (Dependency Inversion Principle) so that
//! specialized analyzers (GPU-accelerated, anisotropic-metric-aware, etc.) can
//! be substituted without changing call sites.
//!
//! [`StandardQualityAnalyzer`] delegates to the existing per-triangle functions
//! in [`super::triangle`] and the [`MeshValidator`] for the base report, then
//! enriches it with histograms and aggregate counts.

use crate::application::quality::histograms::Histogram;
use crate::application::quality::report::FullQualityReport;
use crate::application::quality::triangle;
use crate::application::quality::validation::{MeshValidator, QualityThresholds};
use crate::domain::core::scalar::Point3r;
use crate::domain::mesh::IndexedMesh;

// ── Trait ─────────────────────────────────────────────────────────────────────

/// DIP-compliant interface for computing a [`FullQualityReport`] from an [`IndexedMesh`].
pub trait QualityAnalyzer {
    /// Compute the full quality report for `mesh`.
    fn compute(&self, mesh: &IndexedMesh) -> FullQualityReport;
}

// ── Standard implementation ────────────────────────────────────────────────────

/// Standard quality analyzer: delegates to triangle-level metrics.
///
/// Use [`StandardQualityAnalyzer::default()`] for sensible CFD defaults or
/// construct with custom thresholds and histogram bin count.
pub struct StandardQualityAnalyzer {
    /// Thresholds used for pass/fail decisions and bad-face counts.
    pub thresholds: QualityThresholds,
    /// Number of histogram bins.  Defaults to 20.
    pub n_histogram_bins: usize,
}

impl Default for StandardQualityAnalyzer {
    fn default() -> Self {
        Self {
            thresholds: QualityThresholds::default(),
            n_histogram_bins: 20,
        }
    }
}

impl QualityAnalyzer for StandardQualityAnalyzer {
    fn compute(&self, mesh: &IndexedMesh) -> FullQualityReport {
        let n = mesh.faces.len();
        let mut aspect_values = Vec::with_capacity(n);
        let mut min_angle_values = Vec::with_capacity(n);
        let mut edge_len_values = Vec::with_capacity(n * 3);
        let mut bad_angle_count = 0usize;
        let mut bad_aspect_count = 0usize;

        for face in mesh.faces.iter() {
            let pa = mesh.vertices.position(face.vertices[0]);
            let pb = mesh.vertices.position(face.vertices[1]);
            let pc = mesh.vertices.position(face.vertices[2]);
            let a = Point3r::new(pa.x, pa.y, pa.z);
            let b = Point3r::new(pb.x, pb.y, pb.z);
            let c = Point3r::new(pc.x, pc.y, pc.z);

            let ar = triangle::aspect_ratio(&a, &b, &c);
            let ma = triangle::min_angle(&a, &b, &c);
            let el_ab = (b - a).norm();
            let el_bc = (c - b).norm();
            let el_ca = (a - c).norm();

            aspect_values.push(ar);
            min_angle_values.push(ma);
            edge_len_values.push(el_ab);
            edge_len_values.push(el_bc);
            edge_len_values.push(el_ca);

            if ma < self.thresholds.min_angle {
                bad_angle_count += 1;
            }
            if ar > self.thresholds.max_aspect_ratio {
                bad_aspect_count += 1;
            }
        }

        let base = MeshValidator::with_thresholds(self.thresholds.clone())
            .validate(&mesh.faces, &mesh.vertices);

        let bins = self.n_histogram_bins;
        FullQualityReport {
            base,
            edge_length_histogram: Histogram::compute(&edge_len_values, bins),
            min_angle_histogram: Histogram::compute(&min_angle_values, bins),
            aspect_ratio_histogram: Histogram::compute(&aspect_values, bins),
            bad_angle_count,
            bad_aspect_count,
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, PrimitiveMesh};

    fn unit_cube() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube build")
    }

    #[test]
    fn analyzer_produces_non_empty_report_for_cube() {
        let mesh = unit_cube();
        let analyzer = StandardQualityAnalyzer::default();
        let report = analyzer.compute(&mesh);
        assert!(report.base.total_faces > 0, "cube should have faces");
        assert!(report.aspect_ratio_histogram.is_some());
        assert!(report.min_angle_histogram.is_some());
        assert!(report.edge_length_histogram.is_some());
    }

    #[test]
    fn analyzer_counts_match_face_count() {
        let mesh = unit_cube();
        let analyzer = StandardQualityAnalyzer::default();
        let report = analyzer.compute(&mesh);
        assert_eq!(report.base.total_faces, mesh.faces.len());
    }

    #[test]
    fn cube_has_no_bad_angle_faces() {
        let mesh = unit_cube();
        // Default threshold is 15°; cube right-triangles have min angle ~45°.
        let analyzer = StandardQualityAnalyzer::default();
        let report = analyzer.compute(&mesh);
        assert_eq!(
            report.bad_angle_count, 0,
            "cube triangles have no angles below 15°"
        );
    }

    #[test]
    fn trait_object_dispatch_works() {
        let mesh = unit_cube();
        let analyzer: &dyn QualityAnalyzer = &StandardQualityAnalyzer::default();
        let report = analyzer.compute(&mesh);
        assert!(report.base.total_faces > 0);
    }
}
