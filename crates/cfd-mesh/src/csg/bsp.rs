//! BSP tree for CSG operations.
//!
//! Adapted from csgrs's BSP approach but operating on indexed faces.
//! Splitting plane heuristic: minimize `8 × spans + |front - back|`.

use crate::core::index::{FaceId, RegionId};
use crate::core::scalar::{Real, Point3r, Vector3r};
use crate::geometry::plane::Plane;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::classify::{classify_triangle, PolygonClassification};
use crate::csg::split::split_triangle;

/// Plane-classification tolerance for BSP operations.
///
/// This must be larger than vertex-weld `TOLERANCE` to absorb the floating-
/// point error that accumulates when evaluating `dot(n, p) + w` for points
/// that have passed through multiple arithmetic operations. Using `TOLERANCE`
/// (1e-9 for f64) causes nearly every triangle to be classified `Spanning`,
/// triggering exponential face-splitting and infinite loops.
///
/// 1e-5 mm (~10 nm) matches csgrs's default epsilon and is safe for
/// millimeter-scale geometry.
pub const BSP_PLANE_TOLERANCE: Real = 1e-5;

/// Maximum BSP tree depth to prevent pathological splitting.
const MAX_BSP_DEPTH: usize = 128;

/// A node in the BSP tree.
pub struct BspNode {
    /// Splitting plane for this node.
    plane: Option<Plane>,
    /// Coplanar faces at this node.
    coplanar: Vec<FaceData>,
    /// Front subtree.
    front: Option<Box<BspNode>>,
    /// Back subtree.
    back: Option<Box<BspNode>>,
    /// Depth of this node in the tree.
    depth: usize,
}

impl BspNode {
    /// Create an empty BSP node at depth 0.
    pub fn new() -> Self {
        Self::with_depth(0)
    }

    /// Create an empty BSP node at given depth.
    fn with_depth(depth: usize) -> Self {
        Self {
            plane: None,
            coplanar: Vec::new(),
            front: None,
            back: None,
            depth,
        }
    }

    /// Build a BSP tree from a set of faces.
    pub fn build(faces: &[FaceData], pool: &mut VertexPool) -> Self {
        if faces.is_empty() {
            return Self::new();
        }

        let mut node = Self::new();
        node.add_faces(faces, pool);
        node
    }

    /// Add faces to this BSP node.
    pub fn add_faces(&mut self, faces: &[FaceData], pool: &mut VertexPool) {
        if faces.is_empty() {
            return;
        }

        // Guard: if depth limit is reached, store all remaining faces as
        // coplanar to prevent exponential blowup on degenerate meshes.
        if self.depth >= MAX_BSP_DEPTH {
            self.coplanar.extend_from_slice(faces);
            return;
        }

        // Choose splitting plane if we don't have one yet
        if self.plane.is_none() {
            self.plane = Some(select_splitting_plane(faces, pool));
        }

        let plane = self.plane.as_ref().unwrap();

        let mut front_list = Vec::new();
        let mut back_list = Vec::new();

        for face in faces {
            let (class, classifications) = classify_triangle(face.vertices, pool, plane);
            match class {
                PolygonClassification::Coplanar => {
                    self.coplanar.push(*face);
                }
                PolygonClassification::Front => {
                    front_list.push(*face);
                }
                PolygonClassification::Back => {
                    back_list.push(*face);
                }
                PolygonClassification::Spanning => {
                    let result = split_triangle(face, classifications, plane, pool);
                    front_list.extend(result.front);
                    back_list.extend(result.back);
                }
            }
        }

        let child_depth = self.depth + 1;
        if !front_list.is_empty() {
            let front = self.front.get_or_insert_with(|| Box::new(BspNode::with_depth(child_depth)));
            front.add_faces(&front_list, pool);
        }

        if !back_list.is_empty() {
            let back = self.back.get_or_insert_with(|| Box::new(BspNode::with_depth(child_depth)));
            back.add_faces(&back_list, pool);
        }
    }

    /// Clip a set of faces to the inside of this BSP tree.
    ///
    /// For CSG operations, "inside" means the back half-space of this BSP tree.
    /// Faces entirely in front are discarded; faces in back are kept.
    /// Coplanar faces are kept only if their normal aligns with the plane normal
    /// (i.e., they face the same direction as the "inside" of the plane).
    pub fn clip_faces(&self, faces: &[FaceData], pool: &mut VertexPool) -> Vec<FaceData> {
        let plane = match &self.plane {
            Some(p) => p,
            None => return faces.to_vec(),
        };

        let mut front_list = Vec::new();
        let mut back_list = Vec::new();

        for face in faces {
            let (class, classifications) = classify_triangle(face.vertices, pool, plane);
            match class {
                PolygonClassification::Front => front_list.push(*face),
                PolygonClassification::Back => back_list.push(*face),
                PolygonClassification::Coplanar => {
                    // Critical: Check coplanar face orientation relative to the plane.
                    //
                    // This mirrors csgrs's orient_plane logic for coplanar polygons.
                    // orient_plane takes a point on the polygon's plane and moves it
                    // along the polygon's normal, then checks which side of the BSP
                    // plane it lands on.
                    //
                    // The result is: BSP_normal · polygon_normal > 0 → FRONT
                    //
                    // In csgrs's clip_polygons, coplanar faces go to:
                    // - front_parts if orient_plane returns FRONT (kept when no front child)
                    // - back_parts otherwise (discarded when no back child)
                    //
                    // For a closed solid, faces with normals pointing OUTWARD are on
                    // the exterior surface. When such a face lies on a BSP splitting plane:
                    // - If face normal · plane normal > 0: face is on the "outside" (FRONT)
                    // - If face normal · plane normal <= 0: face is on the "inside" (BACK)
                    let face_normal = compute_face_normal(face, pool);
                    let dot = face_normal.dot(&plane.normal);
                    if dot > 0.0 {
                        // Face normal aligns with BSP plane normal
                        // → this face is on the exterior side of the BSP's solid
                        // → goes to front_list (kept when no front child)
                        front_list.push(*face);
                    } else {
                        // Face normal opposes BSP plane normal (or perpendicular)
                        // → this face is on the interior side of the BSP's solid
                        // → goes to back_list (discarded when no back child)
                        back_list.push(*face);
                    }
                }
                PolygonClassification::Spanning => {
                    let result = split_triangle(face, classifications, plane, pool);
                    front_list.extend(result.front);
                    back_list.extend(result.back);
                }
            }
        }

        let mut result = match &self.front {
            Some(front) => front.clip_faces(&front_list, pool),
            None => front_list,
        };

        if let Some(back) = &self.back {
            result.extend(back.clip_faces(&back_list, pool));
        }
        // If no back subtree, back_list is discarded (clipped away)

        result
    }

    /// Clip this BSP tree to another BSP tree (modifies this tree in place).
    pub fn clip_to(&mut self, other: &BspNode, pool: &mut VertexPool) {
        self.coplanar = other.clip_faces(&self.coplanar, pool);
        if let Some(ref mut front) = self.front {
            front.clip_to(other, pool);
        }
        if let Some(ref mut back) = self.back {
            back.clip_to(other, pool);
        }
    }

    /// Collect all faces in the BSP tree.
    pub fn all_faces(&self) -> Vec<FaceData> {
        let mut result = self.coplanar.clone();
        if let Some(ref front) = self.front {
            result.extend(front.all_faces());
        }
        if let Some(ref back) = self.back {
            result.extend(back.all_faces());
        }
        result
    }

    /// Invert the BSP tree (flip all planes and swap front/back).
    pub fn invert(&mut self) {
        if let Some(ref mut plane) = self.plane {
            *plane = plane.flip();
        }

        for face in &mut self.coplanar {
            face.flip();
        }

        if let Some(ref mut front) = self.front {
            front.invert();
        }
        if let Some(ref mut back) = self.back {
            back.invert();
        }

        std::mem::swap(&mut self.front, &mut self.back);
    }
}

impl Default for BspNode {
    fn default() -> Self {
        Self::new()
    }
}

/// Select the best splitting plane from a set of faces.
///
/// Heuristic: minimize `4 × spanning_count + |front_count - back_count|`.
///
/// The spanning weight is 4 (not the earlier 8) so that axis-aligned
/// centroid-median planes (which produce balanced splits at the cost of some
/// spanning) can beat the degenerate all-back face planes produced by convex
/// tessellated bodies such as spheres.
///
/// In addition to the first 16 face-plane candidates, the function always
/// tries the three axis-aligned planes through the bounding-box centroid of
/// all face centroids.  For convex polyhedra those planes are well-balanced
/// (score ≈ 4×spanning + 0) whereas any face-plane gives score = N (all
/// other faces on one side, 0 spanning).
fn select_splitting_plane(faces: &[FaceData], pool: &VertexPool) -> Plane {
    let candidate_count = faces.len().min(16);

    let mut best_plane = face_plane(&faces[0], pool);
    let mut best_score = i64::MAX;

    // Closure: score a single candidate plane.
    let score_plane = |candidate: &Plane| -> i64 {
        let mut front = 0i64;
        let mut back = 0i64;
        let mut spanning = 0i64;
        for face in faces {
            let (class, _) = classify_triangle(face.vertices, pool, candidate);
            match class {
                PolygonClassification::Front    => front += 1,
                PolygonClassification::Back     => back += 1,
                PolygonClassification::Spanning => spanning += 1,
                PolygonClassification::Coplanar => {}
            }
        }
        4 * spanning + (front - back).abs()
    };

    // Try face-plane candidates.
    for i in 0..candidate_count {
        let candidate = face_plane(&faces[i], pool);
        let score = score_plane(&candidate);
        if score < best_score {
            best_score = score;
            best_plane = candidate;
        }
    }

    // Also try axis-aligned centroid-median planes.
    //
    // For convex polyhedra (e.g. tessellated spheres), every face plane puts
    // all other faces into the back subtree (score = N), while an axis-aligned
    // median plane produces a balanced split (score ≈ 4×spanning).  Without
    // these extra candidates the BSP degenerates into an O(N) linear chain,
    // causing cascade-splitting of near-seam faces and exponential face blowup.
    let mut cx_sum = 0.0_f64;
    let mut cy_sum = 0.0_f64;
    let mut cz_sum = 0.0_f64;
    let n = faces.len() as Real;
    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);
        cx_sum += (a.x + b.x + c.x) as f64;
        cy_sum += (a.y + b.y + c.y) as f64;
        cz_sum += (a.z + b.z + c.z) as f64;
    }
    let mid_x = (cx_sum / (3.0 * n as f64)) as Real;
    let mid_y = (cy_sum / (3.0 * n as f64)) as Real;
    let mid_z = (cz_sum / (3.0 * n as f64)) as Real;
    let axis_candidates = [
        Plane::from_normal_and_point(Vector3r::x(), &Point3r::new(mid_x, 0.0, 0.0)),
        Plane::from_normal_and_point(Vector3r::y(), &Point3r::new(0.0, mid_y, 0.0)),
        Plane::from_normal_and_point(Vector3r::z(), &Point3r::new(0.0, 0.0, mid_z)),
    ];
    for candidate in &axis_candidates {
        let score = score_plane(candidate);
        if score < best_score {
            best_score = score;
            best_plane = *candidate;
        }
    }

    best_plane
}

/// Create a plane from a face's triangle.
fn face_plane(face: &FaceData, pool: &VertexPool) -> Plane {
    let a = pool.position(face.vertices[0]);
    let b = pool.position(face.vertices[1]);
    let c = pool.position(face.vertices[2]);
    Plane::from_three_points(&a, &b, &c)
        .unwrap_or_else(|| Plane::new(Vector3r::z(), 0.0))
}

/// Compute the normal vector of a triangular face.
fn compute_face_normal(face: &FaceData, pool: &VertexPool) -> Vector3r {
    let a = pool.position(face.vertices[0]);
    let b = pool.position(face.vertices[1]);
    let c = pool.position(face.vertices[2]);
    let ab = b - a;
    let ac = c - a;
    let cross = ab.cross(&ac);
    let len = cross.norm();
    if len > 1e-12 {
        cross / len
    } else {
        Vector3r::z() // Degenerate face, return arbitrary normal
    }
}
