pub mod geometry;
pub mod solver;
pub mod validation;

pub use geometry::TrifurcationGeometry3D;
pub use solver::{TrifurcationConfig3D, TrifurcationSolver3D, TrifurcationSolution3D};
pub use validation::{TrifurcationValidator3D, TrifurcationValidationResult3D};

#[cfg(test)]
mod tests;
