//! BSP tree for CSG operations.
//!
//! Adapted from csgrs's BSP approach but operating on indexed faces.
//! Splitting plane heuristic: minimize `8 × spans + |front - back|`.

use crate::core::index::{FaceId, RegionId};
use crate::core::scalar::{Real, Point3r, Vector3r, TOLERANCE};
use crate::geometry::plane::Plane;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::classify::{classify_triangle, PolygonClassification};
use crate::csg::split::split_triangle;

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
}

impl BspNode {
    /// Create an empty BSP node.
    pub fn new() -> Self {
        Self {
            plane: None,
            coplanar: Vec::new(),
            front: None,
            back: None,
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

        if !front_list.is_empty() {
            let front = self.front.get_or_insert_with(|| Box::new(BspNode::new()));
            front.add_faces(&front_list, pool);
        }

        if !back_list.is_empty() {
            let back = self.back.get_or_insert_with(|| Box::new(BspNode::new()));
            back.add_faces(&back_list, pool);
        }
    }

    /// Clip a set of faces to the inside of this BSP tree.
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
                PolygonClassification::Coplanar => front_list.push(*face),
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
/// Heuristic: minimize `8 × spanning_count + |front_count - back_count|`.
fn select_splitting_plane(faces: &[FaceData], pool: &VertexPool) -> Plane {
    let candidate_count = faces.len().min(16); // Limit candidates for performance

    let mut best_plane = face_plane(&faces[0], pool);
    let mut best_score = i64::MAX;

    for i in 0..candidate_count {
        let candidate = face_plane(&faces[i], pool);
        let mut front = 0i64;
        let mut back = 0i64;
        let mut spanning = 0i64;

        for face in faces {
            let (class, _) = classify_triangle(face.vertices, pool, &candidate);
            match class {
                PolygonClassification::Front => front += 1,
                PolygonClassification::Back => back += 1,
                PolygonClassification::Spanning => spanning += 1,
                PolygonClassification::Coplanar => {}
            }
        }

        let score = 8 * spanning + (front - back).abs();
        if score < best_score {
            best_score = score;
            best_plane = candidate;
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
