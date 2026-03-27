//! Selection granularity — controls what level of geometry is selected.

/// The level of geometric detail at which picking operates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SelectionGranularity {
    /// Select entire bodies (scene nodes).
    #[default]
    Body,
    /// Select individual faces.
    Face,
    /// Select individual edges.
    Edge,
    /// Select individual vertices.
    Vertex,
}
