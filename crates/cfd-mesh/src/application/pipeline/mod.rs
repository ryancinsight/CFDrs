//! `NetworkBlueprint → IndexedMesh` pipeline for millifluidic chip designs.
//!
//! Converts a `cfd_schematics::NetworkBlueprint` into watertight boundary-surface
//! `IndexedMesh` objects suitable for CFD simulation and STL/manufacturing output.
//!
//! # Usage
//!
//! ```rust,ignore
//! use cfd_schematics::interface::presets::venturi_chain;
//! use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
//!
//! let bp = venturi_chain("v1", 0.030, 0.004, 0.002);
//! let out = BlueprintMeshPipeline::run(&bp, &PipelineConfig::default()).unwrap();
//! assert!(out.fluid_mesh.is_watertight());
//! ```

pub mod blueprint_mesh;
pub mod constraint;
pub mod region_map;
pub mod topology;
pub mod well_plate;

pub use blueprint_mesh::{BlueprintMeshPipeline, PipelineConfig, PipelineOutput, SegmentCenterline};
pub use constraint::{DiameterConstraintError, InletOutletConstraint, WallClearanceConstraint};
pub use region_map::{BoundaryLabelMap, RegionMap};
pub use topology::{NetworkTopology, TopologyClass};
pub use well_plate::SbsWellPlate96;
