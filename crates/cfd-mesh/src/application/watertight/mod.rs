//! Watertight mesh verification and repair.
//!
//! Critical for CFD: a watertight mesh has no boundary edges, consistent
//! orientation, and no self-intersections.

pub mod check;
pub mod repair;
pub mod seal;

pub use check::{
    assert_watertight, check_watertight,
    WatertightReport,
};
pub use repair::MeshRepair;
