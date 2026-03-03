/// Information about a neighboring process
#[derive(Debug, Clone)]
pub struct NeighborInfo {
    /// Direction of the neighbor
    pub direction: NeighborDirection,
    /// Number of overlapping cells (ghost cell layers)
    pub overlap: usize,
}

/// Direction of neighboring process
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NeighborDirection {
    /// Left neighbor (lower i indices)
    Left,
    /// Right neighbor (higher i indices)
    Right,
    /// Bottom neighbor (lower j indices)
    Bottom,
    /// Top neighbor (higher j indices)
    Top,
    /// Front neighbor (lower k indices)
    Front,
    /// Back neighbor (higher k indices)
    Back,
}
