//! End-to-end design artifact pipeline.
//!
//! Chains: `DesignCandidate` -> `Blueprint` -> `BlueprintMeshPipeline` -> STL / OpenFOAM / SVG / JSON.
//!
//! This module is gated behind the `mesh-export` feature because it depends on
//! `cfd-mesh` for 3-D surface mesh generation and export.
//!
//! # Example
//!
//! ```rust,ignore
//! use cfd_optim::delivery::pipeline::DesignPipeline;
//! use cfd_optim::{SdtOptimizer, OptimMode, SdtWeights};
//!
//! let top5 = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default())
//!     .top_5()
//!     .expect("optimisation failed");
//!
//! let pipeline = DesignPipeline::new("output/designs");
//! let artifacts = pipeline.export_designs(&top5).expect("export failed");
//! for a in &artifacts {
//!     println!("{}: {} verts, {} faces, watertight={}", a.design_id, a.vertex_count, a.face_count, a.watertight);
//! }
//! ```

use std::path::{Path, PathBuf};

use cfd_mesh::application::pipeline::blueprint_mesh::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::domain::core::index::{FaceId, RegionId};
use cfd_mesh::domain::topology::halfedge::PatchType;
use cfd_mesh::infrastructure::io::openfoam::write_openfoam_polymesh;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;

use crate::OptimError;
use crate::RankedDesign;

// ── Output type ──────────────────────────────────────────────────────────────

/// Artifacts produced for a single design.
///
/// All path fields are absolute once the pipeline completes.  The mesh
/// statistics (`vertex_count`, `face_count`, `watertight`) refer to the
/// **fluid** mesh.
pub struct DesignArtifacts {
    /// Candidate identifier (matches `DesignCandidate::id`).
    pub design_id: String,
    /// Path to the binary STL of the fluid-domain mesh.
    pub fluid_stl: PathBuf,
    /// Path to the binary STL of the chip-body mesh, if generated.
    pub chip_stl: Option<PathBuf>,
    /// Path to the OpenFOAM `constant/polyMesh` directory.
    pub openfoam_dir: PathBuf,
    /// Path to the 2-D channel schematic SVG.
    pub schematic_svg: PathBuf,
    /// Path to the interchange-format JSON (channel system).
    pub schematic_json: PathBuf,
    /// Number of deduplicated vertices in the fluid mesh.
    pub vertex_count: usize,
    /// Number of triangular faces in the fluid mesh.
    pub face_count: usize,
    /// Whether the fluid mesh passed the watertight (closed-manifold) check.
    pub watertight: bool,
}

// ── Pipeline ─────────────────────────────────────────────────────────────────

/// Pipeline that converts ranked designs into mesh and schematic artifacts.
///
/// For each [`RankedDesign`] the pipeline:
///
/// 1. Converts the candidate to a [`cfd_schematics::NetworkBlueprint`].
/// 2. Runs [`BlueprintMeshPipeline`] to produce fluid and chip meshes.
/// 3. Writes binary STL files for both meshes.
/// 4. Writes an OpenFOAM `constant/polyMesh` directory with inlet/outlet/wall
///    boundary patches derived from the mesh's boundary labels.
/// 5. Saves a 2-D channel schematic SVG via `cfd-schematics`.
/// 6. Saves the interchange-format JSON for the channel system.
/// 7. Returns a [`DesignArtifacts`] summary with mesh statistics.
pub struct DesignPipeline {
    output_root: PathBuf,
    mesh_config: PipelineConfig,
}

impl DesignPipeline {
    /// Create a new pipeline that writes artifacts under `output_root`.
    ///
    /// Each design gets its own sub-directory named after the candidate ID.
    /// The diameter constraint is skipped by default because millifluidic
    /// designs have channel cross-sections smaller than the 4 mm macro-port
    /// spec (physical tubing adapters are external).
    pub fn new(output_root: impl AsRef<Path>) -> Self {
        let mut config = PipelineConfig::default();
        config.skip_diameter_constraint = true;
        Self {
            output_root: output_root.as_ref().to_path_buf(),
            mesh_config: config,
        }
    }

    /// Create a pipeline with a custom [`PipelineConfig`].
    pub fn with_config(output_root: impl AsRef<Path>, config: PipelineConfig) -> Self {
        Self {
            output_root: output_root.as_ref().to_path_buf(),
            mesh_config: config,
        }
    }

    /// Export all artifacts for a single design.
    ///
    /// # Errors
    ///
    /// Returns [`OptimError::MeshExportFailed`] if any stage of the pipeline
    /// (mesh generation, file I/O, schematic rendering) fails.
    pub fn export_design(&self, design: &RankedDesign) -> Result<DesignArtifacts, OptimError> {
        let candidate = &design.candidate;
        let design_id = &candidate.id;

        // 1. Create output directory
        let design_dir = self.output_root.join(design_id);
        std::fs::create_dir_all(&design_dir).map_err(|e| OptimError::MeshExportFailed {
            scenario_id: design_id.clone(),
            message: format!("failed to create output directory: {e}"),
        })?;

        // 2. Convert candidate to blueprint
        let blueprint = candidate.to_blueprint();

        // 3. Run BlueprintMeshPipeline
        let mut output =
            BlueprintMeshPipeline::run(&blueprint, &self.mesh_config).map_err(|e| {
                OptimError::MeshExportFailed {
                    scenario_id: design_id.clone(),
                    message: format!("mesh pipeline failed: {e}"),
                }
            })?;

        // 4. Write fluid STL (binary)
        let fluid_stl_path = design_dir.join("fluid.stl");
        {
            let mut file =
                std::io::BufWriter::new(std::fs::File::create(&fluid_stl_path).map_err(|e| {
                    OptimError::MeshExportFailed {
                        scenario_id: design_id.clone(),
                        message: format!("failed to create fluid STL file: {e}"),
                    }
                })?);
            write_stl_binary(&mut file, &output.fluid_mesh).map_err(|e| {
                OptimError::MeshExportFailed {
                    scenario_id: design_id.clone(),
                    message: format!("failed to write fluid STL: {e}"),
                }
            })?;
        }

        // 4b. Write chip STL if the pipeline produced one
        let chip_stl_path = if let Some(ref chip_mesh) = output.chip_mesh {
            let path = design_dir.join("chip.stl");
            let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(|e| {
                OptimError::MeshExportFailed {
                    scenario_id: design_id.clone(),
                    message: format!("failed to create chip STL file: {e}"),
                }
            })?);
            write_stl_binary(&mut file, chip_mesh).map_err(|e| OptimError::MeshExportFailed {
                scenario_id: design_id.clone(),
                message: format!("failed to write chip STL: {e}"),
            })?;
            Some(path)
        } else {
            None
        };

        // 5. Write OpenFOAM polyMesh with region remapping from boundary labels
        let openfoam_dir = design_dir.join("constant").join("polyMesh");
        let mut fluid_for_of = output.fluid_mesh.clone();
        reassign_regions_from_labels(&mut fluid_for_of);
        write_openfoam_polymesh(
            &fluid_for_of,
            &openfoam_dir,
            &[
                (RegionId::new(0), "inlet", PatchType::Inlet),
                (RegionId::new(1), "outlet", PatchType::Outlet),
                (RegionId::new(2), "walls", PatchType::Wall),
            ],
        )
        .map_err(|e| OptimError::MeshExportFailed {
            scenario_id: design_id.clone(),
            message: format!("failed to write OpenFOAM polyMesh: {e}"),
        })?;

        // 6. Save schematic SVG
        let schematic_svg_path = design_dir.join("schematic.svg");
        super::export::save_schematic_svg(candidate, &schematic_svg_path).map_err(|e| {
            OptimError::MeshExportFailed {
                scenario_id: design_id.clone(),
                message: format!("failed to save schematic SVG: {e}"),
            }
        })?;

        // 7. Save blueprint-native schematic JSON
        let schematic_json_path = design_dir.join("schematic.json");
        let blueprint = candidate.to_blueprint();
        let json = blueprint
            .to_json_pretty()
            .map_err(|e| OptimError::MeshExportFailed {
                scenario_id: design_id.clone(),
                message: format!("failed to serialize blueprint JSON: {e}"),
            })?;
        std::fs::write(&schematic_json_path, json).map_err(|e| OptimError::MeshExportFailed {
            scenario_id: design_id.clone(),
            message: format!("failed to write schematic JSON file: {e}"),
        })?;

        // 8. Collect mesh statistics
        let vertex_count = output.fluid_mesh.vertex_count();
        let face_count = output.fluid_mesh.face_count();
        let watertight = output.fluid_mesh.is_watertight();

        Ok(DesignArtifacts {
            design_id: design_id.clone(),
            fluid_stl: fluid_stl_path,
            chip_stl: chip_stl_path,
            openfoam_dir,
            schematic_svg: schematic_svg_path,
            schematic_json: schematic_json_path,
            vertex_count,
            face_count,
            watertight,
        })
    }

    /// Export artifacts for multiple designs.
    ///
    /// Processes each design sequentially and collects the results.  Stops on
    /// the first error.
    ///
    /// # Errors
    ///
    /// Returns the first [`OptimError::MeshExportFailed`] encountered.
    pub fn export_designs(
        &self,
        designs: &[RankedDesign],
    ) -> Result<Vec<DesignArtifacts>, OptimError> {
        designs.iter().map(|d| self.export_design(d)).collect()
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Re-assign face `RegionId` values based on the boundary labels produced by
/// the mesh pipeline.
///
/// Mapping:
/// - `"inlet"` -> `RegionId(0)`
/// - `"outlet"` -> `RegionId(1)`
/// - everything else -> `RegionId(2)` (wall)
fn reassign_regions_from_labels(mesh: &mut cfd_mesh::domain::mesh::IndexedMesh) {
    use std::collections::HashMap;

    let region_for: HashMap<FaceId, RegionId> = mesh
        .boundary_labels
        .iter()
        .map(|(&fid, label)| {
            let region = match label.as_str() {
                "inlet" => RegionId::new(0),
                "outlet" => RegionId::new(1),
                _ => RegionId::new(2),
            };
            (fid, region)
        })
        .collect();

    for (fid, face) in mesh.faces.iter_mut_enumerated() {
        face.region = *region_for.get(&fid).unwrap_or(&RegionId::new(2));
    }
}
