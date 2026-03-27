//! Selection target — identifies a specific geometric entity in the scene.

/// A selectable entity in the scene, identified by scene node index and
/// optional sub-element index (face, edge, or vertex within the mesh).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SelectionTarget {
    /// An entire scene node (body-level selection).
    Body(usize),
    /// A specific face within a scene node.
    Face(usize, u32),
    /// A specific edge within a scene node.
    Edge(usize, u32),
    /// A specific vertex within a scene node.
    Vertex(usize, u32),
}

impl SelectionTarget {
    /// The scene node index regardless of granularity.
    #[must_use]
    pub fn node_index(&self) -> usize {
        match self {
            Self::Body(idx) | Self::Face(idx, _) | Self::Edge(idx, _) | Self::Vertex(idx, _) => {
                *idx
            }
        }
    }

    /// The sub-element index, if any (None for body-level selection).
    #[must_use]
    pub fn sub_element(&self) -> Option<u32> {
        match self {
            Self::Body(_) => None,
            Self::Face(_, id) | Self::Edge(_, id) | Self::Vertex(_, id) => Some(*id),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn body_has_no_sub_element() {
        let t = SelectionTarget::Body(0);
        assert_eq!(t.node_index(), 0);
        assert_eq!(t.sub_element(), None);
    }

    #[test]
    fn face_has_sub_element() {
        let t = SelectionTarget::Face(3, 42);
        assert_eq!(t.node_index(), 3);
        assert_eq!(t.sub_element(), Some(42));
    }
}
