//! Symmetric Sparse Voxel Directed Acyclic Graph (SSVDAG)
//!
//! Exposes exact scalar-field partitioning and mathematically-pure CSG operations
//! via structural boolean DAG traversal.

pub mod boolean;
pub mod core;
pub mod rasterize;

pub use core::{DagIndex, SparseVoxelOctree, SvoNode};
