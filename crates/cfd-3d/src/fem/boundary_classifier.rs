//! Axial boundary face classifier for 3D channel FEM meshes.
//!
//! # Algorithm — Z-Axial Boundary Classification for Extruded Meshes
//!
//! For channel geometries extruded along the z-axis, classifies boundary
//! faces into inlet (z ≈ z_min), outlet (z ≈ z_max), and wall (lateral)
//! using either explicit mesh labels or a tolerance derived from the axial
//! mesh resolution.
//!
//! ## Strategy
//!
//! ```text
//! 1. Collect all topological boundary faces (faces belonging to exactly one cell).
//! 2. If ANY face has an explicit label → label-based path:
//!      "inlet"          → inlet_nodes
//!      "outlet" or
//!      "outlet_<N>"     → outlet_nodes (and outlet_nodes_by_label["outlet_<N>"])
//!      other / unlabeled → wall_nodes (if not already in inlet/outlet set)
//! 3. If NO face has a label → z-proximity fallback:
//!      a. Compute global z_min, z_max over all mesh vertices.
//!      b. Tolerance Δz = (z_max − z_min) / n_axial × 0.75  (sub-cell tolerance).
//!      c. Face z-center = mean z of face vertices.
//!      d. |z_center − z_min| ≤ Δz  → inlet
//!         |z_center − z_max| ≤ Δz  → outlet
//!         otherwise                  → wall
//! ```
//!
//! # Theorem — Sufficient Condition for Z-Proximity Classification
//!
//! For a mesh with z ∈ [z_min, z_max] extruded with n_axial cells of uniform
//! size Δz = z_span/n_axial, a face whose centroid satisfies
//! |z_c − z_endpoint| < 0.75 Δz cannot simultaneously be classified as both
//! inlet and outlet (provided z_span > 0 and n_axial ≥ 2), because
//!
//! ```text
//! |z_min − z_max| = z_span ≥ 2 · Δz > 2 · 0.75 · Δz = 1.5 Δz
//! ```
//!
//! so the inlet and outlet tolerance bands are disjoint.
//!
//! ## References
//!
//! - Deardorff, J. W. (1970). "A numerical study of three-dimensional
//!   turbulent channel flow at large Reynolds numbers." *J. Fluid Mech.*
//!   41(2):453–480.
//! - Elman, H., Silvester, D. & Wathen, A. (2014). *Finite Elements and
//!   Fast Iterative Solvers.* Oxford University Press, §3.4.

use cfd_mesh::{domain::core::index::FaceId, domain::core::Scalar, IndexedMesh};
use num_traits::Float;
use std::collections::{HashMap, HashSet};

// ── Named tolerancing constant ────────────────────────────────────────────────

/// Sub-cell tolerance factor for z-proximity classification.
///
/// A face whose z-centroid lies within `AXIAL_TOL_FACTOR × Δz` of z_min or
/// z_max is classified as inlet or outlet respectively.
/// The value 0.75 ensures the tolerance is strictly less than one cell width,
/// keeping the inlet/outlet bands disjoint for n_axial ≥ 2.
const AXIAL_TOL_FACTOR: f64 = 0.75;

// ═══════════════════════════════════════════════════════════════════════════════
// Output type
// ═══════════════════════════════════════════════════════════════════════════════

/// Node index sets produced by [`AxialBoundaryClassifier::classify`].
///
/// All sets use vertex indices in the form `VertexId::as_usize()` so that
/// callers can directly look up `boundary_conditions[&node_idx]`.
#[derive(Debug, Default)]
pub struct BoundaryFaceSets {
    /// Vertex indices whose faces were labelled / classified as "inlet".
    pub inlet_nodes: HashSet<usize>,
    /// Vertex indices whose faces were labelled / classified as "outlet" or "outlet_*".
    ///
    /// This is the union of all outlet labels; use [`outlet_nodes_by_label`] when
    /// per-outlet pressure assignment is needed (e.g., trifurcation).
    pub outlet_nodes: HashSet<usize>,
    /// Per-label outlet node sets.
    ///
    /// For single-outlet meshes this contains at most one entry.
    /// For multi-outlet meshes (bifurcation: "outlet_0"/"outlet_1",
    /// trifurcation: "outlet_0"/"outlet_1"/"outlet_2") each label maps
    /// to its own vertex set so solvers can assign per-outlet pressures.
    pub outlet_nodes_by_label: HashMap<String, HashSet<usize>>,
    /// Vertex indices on wall faces.
    ///
    /// This set MAY overlap with `inlet_nodes` and `outlet_nodes` at *rim nodes* —
    /// vertices shared between an inlet/outlet face and a lateral wall face.
    /// Callers that assign wall BCs should use `.or_insert` (or an explicit
    /// `contains` check) so that inlet/outlet BCs take priority over wall at rim nodes.
    pub wall_nodes: HashSet<usize>,
    /// Union of all boundary vertex indices (inlet ∪ outlet ∪ wall).
    pub boundary_vertices: HashSet<usize>,
}

// ═══════════════════════════════════════════════════════════════════════════════
// Classifier
// ═══════════════════════════════════════════════════════════════════════════════

/// Classifies boundary faces of extruded-channel `IndexedMesh` objects into
/// inlet, outlet, and wall node sets.
///
/// See [module-level documentation](self) for the full algorithm.
///
/// # Usage
///
/// ```rust,ignore
/// let classifier = AxialBoundaryClassifier::new(&mesh, config.resolution.0);
/// let face_sets  = classifier.classify();
/// // Apply inlet velocity to face_sets.inlet_nodes …
/// ```
pub struct AxialBoundaryClassifier<'a, T: Scalar> {
    mesh: &'a IndexedMesh<T>,
    /// Number of axial (z-direction) cells — sets the z-proximity tolerance.
    resolution_axial: usize,
}

impl<'a, T: Scalar> AxialBoundaryClassifier<'a, T> {
    /// Create a classifier for `mesh`.
    ///
    /// `resolution_axial` is the number of cells along the extrusion axis;
    /// used to compute the z-proximity tolerance for any unlabeled faces.
    #[must_use]
    pub fn new(mesh: &'a IndexedMesh<T>, resolution_axial: usize) -> Self {
        Self {
            mesh,
            resolution_axial,
        }
    }

    /// Classify all boundary faces and return the partitioned node sets.
    ///
    /// # Algorithm — Unified Hybrid Classification
    ///
    /// For each boundary face, the label takes priority; for unlabeled faces
    /// z-proximity to the global z-endpoints is used as a fallback.  This
    /// handles meshes where:
    ///  - all faces are labeled (bifurcation, serpentine, trifurcation),
    ///  - no faces are labeled (pure z-proximity),
    ///  - or lateral/wall faces are labeled while inlet/outlet are unlabeled
    ///    (venturi meshes from VenturiMeshBuilder).
    ///
    /// ```text
    /// Pass 1 — inlet and outlet:
    ///   For each boundary face f:
    ///     label = mesh.boundary_label(f)
    ///     if label == "inlet"          → inlet_nodes   ∪= vertices(f)
    ///     if label starts "outlet"     → outlet_nodes  ∪= vertices(f)
    ///     if label is None:
    ///       z_center = mean z of face vertices
    ///       if |z_center − z_min| ≤ z_tol  → inlet_nodes  ∪= vertices(f)
    ///       if |z_center − z_max| ≤ z_tol  → outlet_nodes ∪= vertices(f)
    ///
    /// Pass 2 — wall (all non-inlet/non-outlet face vertices):
    ///   For each boundary face f not labeled "inlet" or "outlet*":
    ///     wall_nodes ∪= vertices(f)
    ///   (Rim nodes shared with inlet/outlet faces appear in BOTH wall_nodes
    ///    AND inlet/outlet_nodes; callers use .or_insert for wall BCs so that
    ///    inlet/outlet BCs take priority at rim nodes.)
    /// ```
    #[must_use]
    pub fn classify(&self) -> BoundaryFaceSets {
        let mut sets = BoundaryFaceSets::default();
        let boundary_faces = self.mesh.boundary_faces();

        // Pre-compute z-bounds (needed for unlabeled face fallback).
        let (z_min, z_max, z_tol) = self.compute_z_bounds();

        self.classify_unified(&boundary_faces, &mut sets, z_min, z_max, z_tol);
        sets
    }

    // ── Z-bounds helper ───────────────────────────────────────────────────────

    fn compute_z_bounds(&self) -> (T, T, T) {
        let mut z_min = <T as nalgebra::RealField>::max_value().unwrap_or_else(T::one);
        let mut z_max = -<T as nalgebra::RealField>::max_value().unwrap_or_else(T::one);
        for (_, v) in self.mesh.vertices.iter() {
            let z = v.position.z;
            if z < z_min {
                z_min = z;
            }
            if z > z_max {
                z_max = z;
            }
        }
        let z_span = z_max - z_min;
        let n_axial = <T as Scalar>::from_f64(self.resolution_axial.max(1) as f64);
        let z_tol = z_span / n_axial * <T as Scalar>::from_f64(AXIAL_TOL_FACTOR);
        (z_min, z_max, z_tol)
    }

    // ── Unified hybrid classification (label + z-proximity fallback) ──────────

    fn classify_unified(
        &self,
        boundary_faces: &[FaceId],
        sets: &mut BoundaryFaceSets,
        z_min: T,
        z_max: T,
        z_tol: T,
    ) {
        // Pass 1: inlet and outlet (labeled OR z-proximity for unlabeled faces).
        for &f_id in boundary_faces {
            let face = self.mesh.faces.get(f_id);
            let label = self.mesh.boundary_label(f_id);

            for &v_idx in &face.vertices {
                sets.boundary_vertices.insert(v_idx.as_usize());
            }

            match label {
                Some("inlet") => {
                    for &v_idx in &face.vertices {
                        sets.inlet_nodes.insert(v_idx.as_usize());
                    }
                }
                Some(l) if l == "outlet" || l.starts_with("outlet_") => {
                    let per_label = sets.outlet_nodes_by_label.entry(l.to_string()).or_default();
                    for &v_idx in &face.vertices {
                        let id = v_idx.as_usize();
                        sets.outlet_nodes.insert(id);
                        per_label.insert(id);
                    }
                }
                None => {
                    // Unlabeled face: classify by z-proximity to endpoints.
                    let z_center = self.face_z_center(face);
                    if Float::abs(z_center - z_min) <= z_tol {
                        for &v_idx in &face.vertices {
                            sets.inlet_nodes.insert(v_idx.as_usize());
                        }
                    } else if Float::abs(z_center - z_max) <= z_tol {
                        let per_label = sets
                            .outlet_nodes_by_label
                            .entry("outlet".to_string())
                            .or_default();
                        for &v_idx in &face.vertices {
                            let id = v_idx.as_usize();
                            sets.outlet_nodes.insert(id);
                            per_label.insert(id);
                        }
                    }
                    // else: unlabeled AND not near endpoints → handled in pass 2
                }
                Some(_) => {} // e.g., "wall" → handled in pass 2
            }
        }

        // Pass 2: wall — all vertices from non-inlet / non-outlet faces.
        //
        // Rim nodes (vertices shared between an inlet/outlet face and a lateral
        // wall face) are intentionally included in `wall_nodes` even though they
        // are also in `inlet_nodes` / `outlet_nodes`.  This allows callers to
        // detect rim nodes via set intersection and assign appropriate BCs
        // (e.g., no-slip Dirichlet at the inlet rim in the Venturi solver).
        // Callers that do not need rim detection should use `.or_insert` when
        // applying wall BCs so the higher-priority inlet/outlet BC is preserved.
        for &f_id in boundary_faces {
            let face = self.mesh.faces.get(f_id);
            let label = self.mesh.boundary_label(f_id);

            // Skip faces that are explicitly labeled as inlet or outlet — those
            // vertices are already in the inlet/outlet sets with full priority.
            let is_inlet_face = matches!(label, Some("inlet"));
            let is_outlet_face =
                matches!(label, Some(l) if l == "outlet" || l.starts_with("outlet_"));
            if is_inlet_face || is_outlet_face {
                continue;
            }

            for &v_idx in &face.vertices {
                sets.wall_nodes.insert(v_idx.as_usize());
            }
        }
    }

    // ── Per-face z-centroid helper ────────────────────────────────────────────

    fn face_z_center(&self, face: &cfd_mesh::infrastructure::storage::face_store::FaceData) -> T {
        let mut z_sum = T::zero();
        let mut count = 0usize;
        for &v_idx in &face.vertices {
            z_sum = z_sum + self.mesh.vertices.get(v_idx).position.z;
            count += 1;
        }
        if count > 0 {
            z_sum / <T as Scalar>::from_f64(count as f64)
        } else {
            T::zero()
        }
    }
}
