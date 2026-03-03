/// Domain decomposition strategy
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompositionStrategy {
    /// Simple 1D decomposition along x-direction
    Simple1D,
    /// 2D Cartesian decomposition (x and y directions)
    Cartesian2D,
    /// Recursive coordinate bisection (for load balancing)
    RecursiveBisection,
    /// METIS-based graph partitioning (if available)
    Metis,
}
