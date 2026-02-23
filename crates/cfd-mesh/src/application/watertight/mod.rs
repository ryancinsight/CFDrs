//! Watertight mesh verification and repair.
//!
//! Critical for CFD: a watertight mesh has no boundary edges, consistent
//! orientation, and no self-intersections.

pub mod check;
pub mod repair;
pub mod seal;

pub use check::{WatertightReport, check_watertight, assert_watertight,
                check_halfedge, assert_halfedge_watertight};
pub use repair::MeshRepair;
