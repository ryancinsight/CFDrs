//! CFD volume-cell quality metrics.
//!
//! Two metrics that OpenFOAM reports for every internal face:
//!
//! | Metric | Formula | Good | Acceptable | Invalid |
//! |--------|---------|------|------------|---------|
//! | Non-orthogonality | ∠(face normal, d) | < 70° | < 85° | ≥ 90° |
//! | Skewness | \|fc - Pi\| / \|fc - C_owner\| | < 0.5 | < 0.85 | ≥ 1.0 |
//!
//! where **d** is the owner→neighbour centroid vector, **fc** is the face
//! centre, and **Pi** is the point where **d** intersects the face plane.

use nalgebra::Point3;

use crate::application::quality::metrics::QualityMetric;
use crate::domain::core::index::{CellId, FaceId};
use crate::domain::core::scalar::Real;
use crate::domain::mesh::IndexedMesh;

// ── Per-face computations ─────────────────────────────────────────────────────

/// Non-orthogonality in degrees: angle between the face normal and the
/// owner→neighbour centroid vector **d**.
///
/// Returns `0.0` when `mesh` has no cells.
pub fn face_non_orthogonality(
    face: FaceId,
    owner: CellId,
    neighbour: CellId,
    mesh: &IndexedMesh,
) -> Real {
    let face_data = mesh.faces.get(face);
    let [va, vb, vc] = face_data.vertices;
    let a = mesh.vertices.position(va);
    let b = mesh.vertices.position(vb);
    let c = mesh.vertices.position(vc);

    let n = (b - a).cross(&(c - a));
    if n.norm_squared() < 1e-30 {
        return 0.0;
    }
    let n = n.normalize();

    let c_owner = cell_centroid(owner, mesh);
    let c_neigh  = cell_centroid(neighbour, mesh);
    let d = c_neigh - c_owner;
    if d.norm_squared() < 1e-30 {
        return 0.0;
    }
    let d = d.normalize();

    let cos_theta = n.dot(&d).abs().min(1.0);
    cos_theta.acos().to_degrees()
}

/// Skewness: ratio of (face-centre deviation from the owner→neighbour
/// intersection point) to (distance from face centre to owner centroid).
///
/// Returns `0.0` when `mesh` has no cells or the geometry is degenerate.
pub fn face_skewness(
    face: FaceId,
    owner: CellId,
    neighbour: CellId,
    mesh: &IndexedMesh,
) -> Real {
    let face_data = mesh.faces.get(face);
    let [va, vb, vc] = face_data.vertices;
    let a = mesh.vertices.position(va);
    let b = mesh.vertices.position(vb);
    let c = mesh.vertices.position(vc);

    // Face centre = centroid of the triangle.
    let fc = Point3::from((a.coords + b.coords + c.coords) / 3.0);

    // Face normal (unnormalised; used for plane intersection).
    let n = (b - a).cross(&(c - a));
    if n.norm_squared() < 1e-30 {
        return 0.0;
    }

    let c_owner = cell_centroid(owner, mesh);
    let c_neigh  = cell_centroid(neighbour, mesh);

    // Parametric intersection of the d-line with the face plane:
    //   P(t) = c_owner + t * d,  and  n · (P(t) - fc) = 0
    let d = c_neigh - c_owner;
    let denom = n.dot(&d);
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    let t = n.dot(&(fc - c_owner)) / denom;
    let p_i = c_owner + d * t; // intersection point on face plane

    let deviation = (fc - p_i).norm();
    let ref_dist  = (fc - c_owner).norm();
    if ref_dist < 1e-30 {
        return 0.0;
    }
    deviation / ref_dist
}

// ── Report ────────────────────────────────────────────────────────────────────

/// Summary of volume-cell quality for an entire [`IndexedMesh`].
#[derive(Clone, Debug)]
pub struct CellQualityReport {
    /// Non-orthogonality statistics (degrees) over all internal faces.
    pub non_orthogonality: QualityMetric,
    /// Skewness statistics over all internal faces.
    pub skewness: QualityMetric,
    /// Number of internal faces with non-orthogonality > 70°.
    pub high_non_orthogonality_count: usize,
    /// Number of internal faces with skewness > 0.85.
    pub high_skewness_count: usize,
    /// Total internal faces evaluated.
    pub internal_face_count: usize,
}

/// Compute cell quality metrics for all internal faces.
///
/// Returns `None` when the mesh has no volumetric cells or no internal faces.
pub fn cell_quality_report(mesh: &IndexedMesh) -> Option<CellQualityReport> {
    if mesh.cell_count() == 0 {
        return None;
    }

    // Build face → (owner, optional neighbour) map.
    let mut face_owner:    std::collections::HashMap<FaceId, CellId> = std::collections::HashMap::new();
    let mut face_neighbour: std::collections::HashMap<FaceId, CellId> = std::collections::HashMap::new();

    for (cell_id, cell) in mesh.cells_iter_enumerated() {
        for &fi in &cell.faces {
            if let std::collections::hash_map::Entry::Vacant(e) = face_owner.entry(fi) {
                e.insert(cell_id);
            } else {
                face_neighbour.insert(fi, cell_id);
            }
        }
    }

    let mut non_orth_vals: Vec<Real> = Vec::new();
    let mut skew_vals:     Vec<Real> = Vec::new();

    for (&fi, &owner) in &face_owner {
        let Some(&neighbour) = face_neighbour.get(&fi) else { continue };
        non_orth_vals.push(face_non_orthogonality(fi, owner, neighbour, mesh));
        skew_vals.push(face_skewness(fi, owner, neighbour, mesh));
    }

    let non_orthogonality = QualityMetric::from_values(&non_orth_vals)?;
    let skewness          = QualityMetric::from_values(&skew_vals)?;
    let high_no  = non_orth_vals.iter().filter(|&&v| v > 70.0).count();
    let high_sk  = skew_vals.iter().filter(|&&v| v > 0.85).count();

    Some(CellQualityReport {
        non_orthogonality,
        skewness,
        high_non_orthogonality_count: high_no,
        high_skewness_count: high_sk,
        internal_face_count: non_orth_vals.len(),
    })
}

// ── Helper ────────────────────────────────────────────────────────────────────

/// Centroid of a cell: arithmetic mean of its vertices.
pub fn cell_centroid(cell_id: CellId, mesh: &IndexedMesh) -> Point3<Real> {
    let cell = mesh.cell(cell_id);

    // Prefer vertex_ids if populated; fall back to face-vertex union.
    if !cell.vertex_ids.is_empty() {
        let sum: nalgebra::Vector3<Real> = cell.vertex_ids.iter()
            .map(|&vi| mesh.vertices.position(vi).coords)
            .sum();
        return Point3::from(sum / cell.vertex_ids.len() as Real);
    }

    let mut sum = nalgebra::Vector3::<Real>::zeros();
    let mut count = 0usize;
    let mut seen: Vec<_> = Vec::new();
    for &fi in &cell.faces {
        for &vi in &mesh.faces.get(fi).vertices {
            if !seen.contains(&vi) {
                seen.push(vi);
                sum += mesh.vertices.position(vi).coords;
                count += 1;
            }
        }
    }
    if count == 0 {
        Point3::origin()
    } else {
        Point3::from(sum / count as Real)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::grid::StructuredGridBuilder;

    #[test]
    fn cell_quality_report_none_for_surface_mesh() {
        let mesh = IndexedMesh::<f64>::new();
        assert!(cell_quality_report(&mesh).is_none());
    }

    #[test]
    fn regular_grid_non_orthogonality_is_finite() {
        // The 5-tet-per-hex decomposition creates diagonal faces; non-orthogonality
        // can be up to 90° on such faces.  This test verifies the metric is
        // computed correctly (finite, in range [0°, 90°]).
        let mesh = StructuredGridBuilder::new(2, 2, 2).build().unwrap();
        let report = cell_quality_report(&mesh).expect("volume mesh should have report");
        assert!(report.non_orthogonality.max.is_finite());
        assert!(report.non_orthogonality.min >= 0.0);
        assert!(report.non_orthogonality.max <= 90.0 + 1e-9);
        assert!(report.skewness.min >= 0.0);
    }

    #[test]
    fn cell_quality_report_has_internal_faces() {
        let mesh = StructuredGridBuilder::new(2, 2, 2).build().unwrap();
        let report = cell_quality_report(&mesh).unwrap();
        assert!(
            report.internal_face_count > 0,
            "2×2×2 grid should have internal faces"
        );
    }
}
