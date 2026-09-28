//! Two-way and three-way branch junction models with full conservation equations.
//!
//! This module implements the junction flow models governing pressure and flow distribution
//! at branching points in vascular networks. All equations are derived from first principles
//! with references to literature validation.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`two_way_junction`] | `TwoWayBranchJunction<T>` — bifurcation solver |
//! | [`two_way_solution`] | `TwoWayBranchSolution<T>` — solution type |
//! | [`three_way_junction`] | `ThreeWayBranchJunction<T>` + `ThreeWayBranchSolution<T>` |

mod pressure_balance;
pub mod three_way_junction;
pub mod two_way_junction;
pub mod two_way_solution;

pub use three_way_junction::{ThreeWayBranchJunction, ThreeWayBranchSolution};
pub use two_way_junction::TwoWayBranchJunction;
pub use two_way_solution::TwoWayBranchSolution;

#[cfg(test)]
mod tests;
