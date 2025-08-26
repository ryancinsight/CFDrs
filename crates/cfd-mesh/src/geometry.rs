//! Geometric operations for meshes.

/// Geometry operations
pub struct Geometry;
impl Default for Geometry {
    fn default() -> Self {
        Self::new()
    }
}
impl Geometry {
    /// Create a new geometry operations instance
    #[must_use]
    pub fn new() -> Self {
        Self
