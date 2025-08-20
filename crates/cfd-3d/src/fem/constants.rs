//! Constants for FEM fluid dynamics

/// Small value for numerical stability
pub const EPSILON: f64 = 1e-10;

/// Default stabilization parameter
pub const DEFAULT_STABILIZATION: f64 = 0.1;

/// Number of velocity components in 3D
pub const VELOCITY_COMPONENTS: usize = 3;

/// Number of nodes in linear tetrahedron
pub const TET4_NODES: usize = 4;

/// Volume factor for tetrahedron (1/6)
pub const TET_VOLUME_FACTOR: f64 = 6.0;

/// Gauss quadrature points for tetrahedron
pub const GAUSS_POINTS_TET: usize = 4;

/// Shape function tolerance
pub const SHAPE_FUNCTION_TOL: f64 = 1e-12;