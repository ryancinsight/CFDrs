//! OpenFOAM polyMesh writer.
//!
//! Converts an [`IndexedMesh`] or [`HalfEdgeMesh`] into an OpenFOAM
//! `constant/polyMesh` directory that can be used directly with
//! `simpleFoam`, `icoFoam`, or `snappyHexMesh`.
//!
//! ## Output structure
//!
//! ```text
//! <dir>/
//!   points      — vertex coordinates  (vectorField)
//!   faces       — per-face vertex lists (faceList)
//!   owner       — owning-cell index per face (labelList)
//!   neighbour   — internal-face neighbour cells (empty for a surface mesh)
//!   boundary    — named boundary patch definitions (polyBoundaryMesh)
//! ```
//!
//! For a **surface mesh** (closed triangulated shell), every face is a
//! boundary face so the `neighbour` list is empty and the `owner` list
//! simply stores `[0; n_faces]` (each face constitutes its own "cell" in
//! the degenerate 2D sense).  This format is accepted directly by
//! `snappyHexMesh` as the surface geometry input.
//!
//! ## PatchType → OpenFOAM mapping
//!
//! | [`PatchType`] variant | OF `type` string | typical BC |
//! |----------------------|-----------------|------------|
//! | `Inlet`              | `patch`         | fixedValue velocity |
//! | `Outlet`             | `patch`         | fixedValue pressure |
//! | `Wall`               | `wall`          | noSlip |
//! | `Symmetry`           | `symmetry`      | symmetry |
//! | `Periodic`           | `cyclicAMI`     | matchedPair |
//! | `Channel`            | `patch`         | millifluidic channel |
//! | `Custom(_)`          | `patch`         | user-defined |
//!
//! ## Example
//!
//! ```rust,ignore
//! use std::path::Path;
//! use cfd_mesh::{MeshBuilder, io::openfoam::write_openfoam_polymesh};
//! use cfd_mesh::domain::core::index::RegionId;
//! use cfd_mesh::domain::topology::halfedge::PatchType;
//!
//! let mesh = MeshBuilder::new().build(); // populate mesh first
//! write_openfoam_polymesh(
//!     &mesh,
//!     Path::new("constant/polyMesh"),
//!     &[
//!         (RegionId::from_usize(0), "inlet",  PatchType::Inlet),
//!         (RegionId::from_usize(1), "outlet", PatchType::Outlet),
//!         (RegionId::from_usize(2), "walls",  PatchType::Wall),
//!     ],
//! ).expect("OpenFOAM export");
//! ```

use std::io::Write;
use std::path::Path;

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::RegionId;
use crate::domain::mesh::IndexedMesh;
use crate::domain::topology::halfedge::PatchType;

// ── FoamFile header helper ────────────────────────────────────────────────────

/// Write the standard FoamFile header block.
fn write_foam_header(
    w: &mut impl Write,
    class: &str,
    location: &str,
    object: &str,
) -> std::io::Result<()> {
    writeln!(
        w,
        "/*--------------------------------*- C++ -*----------------------------------*\\"
    )?;
    writeln!(
        w,
        "| =========                 |                                                 |"
    )?;
    writeln!(
        w,
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"
    )?;
    writeln!(
        w,
        "|  \\\\    /   O peration     | Version:  v2306                                 |"
    )?;
    writeln!(
        w,
        "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |"
    )?;
    writeln!(
        w,
        "|    \\\\/     M anipulation  |                                                 |"
    )?;
    writeln!(
        w,
        "\\*---------------------------------------------------------------------------*/"
    )?;
    writeln!(w, "FoamFile")?;
    writeln!(w, "{{")?;
    writeln!(w, "    version     2.0;")?;
    writeln!(w, "    format      ascii;")?;
    writeln!(w, "    class       {class};")?;
    writeln!(w, "    location    \"{location}\";")?;
    writeln!(w, "    object      {object};")?;
    writeln!(w, "}}")?;
    writeln!(
        w,
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    )?;
    writeln!(w)?;
    Ok(())
}

// ── Patch descriptor ──────────────────────────────────────────────────────────

/// A single OpenFOAM boundary patch specification.
///
/// Groups a contiguous range of faces under a named boundary with an
/// associated [`PatchType`].  One `PatchSpec` is written per entry in the
/// `boundary` file.
#[derive(Debug, Clone)]
pub struct PatchSpec {
    /// Human-readable patch name used in boundary-condition dictionaries.
    pub name: String,
    /// Physical type of the boundary.
    pub patch_type: PatchType,
    /// First face index belonging to this patch (0-based, within the
    /// re-sorted face array written to the `faces` file).
    pub start_face: usize,
    /// Number of faces in this patch.
    pub n_faces: usize,
}

// ── Public API — IndexedMesh ──────────────────────────────────────────────────

/// Write an [`IndexedMesh`] as an OpenFOAM `constant/polyMesh` directory.
///
/// The files `points`, `faces`, `owner`, `neighbour`, and `boundary` are
/// written directly into `dir`.  The directory is created automatically if it
/// does not exist.
///
/// Faces are re-sorted so that all faces belonging to the same patch are
/// contiguous in the output `faces` / `owner` lists; this is required by
/// the OpenFOAM format specification.
///
/// # Arguments
///
/// - `mesh` — surface mesh to export.
/// - `dir` — output directory (e.g. `Path::new("constant/polyMesh")`).
/// - `patches` — mapping from `RegionId` to `(name, PatchType)`.  Faces
///   whose `RegionId` does not appear in this list are gathered into a
///   synthetic `"defaultFaces"` wall patch.
///
/// # Errors
///
/// Returns [`MeshError::Io`] if any file operation fails, or
/// [`MeshError::Other`] if the mesh has zero vertices.
///
/// # Theorem (face-sort correctness)
///
/// After sorting, the `boundary` patch ranges `[startFace, startFace+nFaces)`
/// are a partition of `[0, nFaces)`.  Every face appears in exactly one
/// patch.  **Proof**: the sort partitions faces by their `RegionId`; unknown
/// IDs fall into the default-faces bucket; the buckets are non-overlapping
/// and union to the full face set. ∎
pub fn write_openfoam_polymesh(
    mesh: &IndexedMesh,
    dir: &Path,
    patches: &[(RegionId, &str, PatchType)],
) -> MeshResult<()> {
    if mesh.vertex_count() == 0 {
        return Err(MeshError::Other(
            "cannot write empty mesh to OpenFOAM format".into(),
        ));
    }

    std::fs::create_dir_all(dir).map_err(MeshError::Io)?;

    // ── Build sorted face list ─────────────────────────────────────────────
    // Build a map: RegionId → (name_index, name, patch_type)
    let patch_map: Vec<(RegionId, String, PatchType)> = patches
        .iter()
        .map(|(rid, name, pt)| (*rid, (*name).to_owned(), pt.clone()))
        .collect();

    // Assign each face to a patch bucket
    let mut buckets: Vec<Vec<[u32; 3]>> = vec![Vec::new(); patch_map.len() + 1]; // +1 for default
    let default_idx = patch_map.len();

    for (_, face) in mesh.faces.iter_enumerated() {
        let bucket = patch_map
            .iter()
            .position(|(rid, _, _)| *rid == face.region)
            .unwrap_or(default_idx);
        buckets[bucket].push([
            face.vertices[0].raw(),
            face.vertices[1].raw(),
            face.vertices[2].raw(),
        ]);
    }

    // Build PatchSpec list (skip empty patches unless it is the default)
    let mut specs: Vec<PatchSpec> = Vec::new();
    let mut cursor = 0usize;
    for (i, (_, name, pt)) in patch_map.iter().enumerate() {
        let n = buckets[i].len();
        if n > 0 {
            specs.push(PatchSpec {
                name: name.clone(),
                patch_type: pt.clone(),
                start_face: cursor,
                n_faces: n,
            });
            cursor += n;
        }
    }
    if !buckets[default_idx].is_empty() {
        let n = buckets[default_idx].len();
        specs.push(PatchSpec {
            name: "defaultFaces".into(),
            patch_type: PatchType::Wall,
            start_face: cursor,
            n_faces: n,
        });
        // No cursor advance needed; this is the last bucket.
    }

    // Flatten sorted face list
    let sorted_faces: Vec<[u32; 3]> = buckets.iter().flatten().copied().collect();

    // ── Write all five polyMesh files ─────────────────────────────────────
    write_points_file(dir, mesh.vertices.positions().copied())?;
    write_faces_file(dir, &sorted_faces)?;
    write_owner_file(dir, sorted_faces.len())?;
    write_neighbour_file(dir)?;
    write_boundary_file(dir, &specs)?;

    Ok(())
}

// ── Shared file writers ───────────────────────────────────────────────────────

fn write_points_file(
    dir: &Path,
    positions: impl Iterator<Item = nalgebra::Point3<crate::domain::core::scalar::Real>>,
) -> MeshResult<()> {
    let path = dir.join("points");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(MeshError::Io)?);
    write_foam_header(&mut file, "vectorField", "constant/polyMesh", "points")
        .map_err(MeshError::Io)?;
    let pts: Vec<_> = positions.collect();
    writeln!(file, "{}", pts.len()).map_err(MeshError::Io)?;
    writeln!(file, "(").map_err(MeshError::Io)?;
    for p in &pts {
        writeln!(file, "    ({:.15e} {:.15e} {:.15e})", p.x, p.y, p.z).map_err(MeshError::Io)?;
    }
    writeln!(file, ")").map_err(MeshError::Io)?;
    writeln!(
        file,
        "\n// ************************************************************************* //"
    )
    .map_err(MeshError::Io)?;
    Ok(())
}

fn write_faces_file(dir: &Path, faces: &[[u32; 3]]) -> MeshResult<()> {
    let path = dir.join("faces");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(MeshError::Io)?);
    write_foam_header(&mut file, "faceList", "constant/polyMesh", "faces")
        .map_err(MeshError::Io)?;
    writeln!(file, "{}", faces.len()).map_err(MeshError::Io)?;
    writeln!(file, "(").map_err(MeshError::Io)?;
    for [v0, v1, v2] in faces {
        writeln!(file, "    3({v0} {v1} {v2})").map_err(MeshError::Io)?;
    }
    writeln!(file, ")").map_err(MeshError::Io)?;
    writeln!(
        file,
        "\n// ************************************************************************* //"
    )
    .map_err(MeshError::Io)?;
    Ok(())
}

fn write_owner_file(dir: &Path, n_faces: usize) -> MeshResult<()> {
    let path = dir.join("owner");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(MeshError::Io)?);
    write_foam_header(&mut file, "labelList", "constant/polyMesh", "owner")
        .map_err(MeshError::Io)?;
    writeln!(file, "{n_faces}").map_err(MeshError::Io)?;
    writeln!(file, "(").map_err(MeshError::Io)?;
    for i in 0..n_faces {
        writeln!(file, "    {i}").map_err(MeshError::Io)?;
    }
    writeln!(file, ")").map_err(MeshError::Io)?;
    writeln!(
        file,
        "\n// ************************************************************************* //"
    )
    .map_err(MeshError::Io)?;
    Ok(())
}

fn write_neighbour_file(dir: &Path) -> MeshResult<()> {
    let path = dir.join("neighbour");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(MeshError::Io)?);
    write_foam_header(&mut file, "labelList", "constant/polyMesh", "neighbour")
        .map_err(MeshError::Io)?;
    writeln!(file, "0").map_err(MeshError::Io)?;
    writeln!(file, "(").map_err(MeshError::Io)?;
    writeln!(file, ")").map_err(MeshError::Io)?;
    writeln!(
        file,
        "\n// ************************************************************************* //"
    )
    .map_err(MeshError::Io)?;
    Ok(())
}

fn write_boundary_file(dir: &Path, specs: &[PatchSpec]) -> MeshResult<()> {
    let path = dir.join("boundary");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path).map_err(MeshError::Io)?);
    write_foam_header(
        &mut file,
        "polyBoundaryMesh",
        "constant/polyMesh",
        "boundary",
    )
    .map_err(MeshError::Io)?;
    writeln!(file, "{}", specs.len()).map_err(MeshError::Io)?;
    writeln!(file, "(").map_err(MeshError::Io)?;
    for spec in specs {
        let of_type = spec.patch_type.openfoam_type();
        writeln!(file, "    {}", spec.name).map_err(MeshError::Io)?;
        writeln!(file, "    {{").map_err(MeshError::Io)?;
        writeln!(file, "        type            {of_type};").map_err(MeshError::Io)?;
        match &spec.patch_type {
            PatchType::Inlet => {
                writeln!(file, "        physicalType    inlet;").map_err(MeshError::Io)?;
            }
            PatchType::Outlet => {
                writeln!(file, "        physicalType    outlet;").map_err(MeshError::Io)?;
            }
            _ => {}
        }
        writeln!(file, "        nFaces          {};", spec.n_faces).map_err(MeshError::Io)?;
        writeln!(file, "        startFace       {};", spec.start_face).map_err(MeshError::Io)?;
        writeln!(file, "    }}").map_err(MeshError::Io)?;
    }
    writeln!(file, ")").map_err(MeshError::Io)?;
    writeln!(
        file,
        "\n// ************************************************************************* //"
    )
    .map_err(MeshError::Io)?;
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::mesh::MeshBuilder;

    /// Helper: build a tiny tetrahedron IndexedMesh.
    fn tet_mesh() -> IndexedMesh {
        let mut b = MeshBuilder::new();
        let a = b.vertex(Point3r::new(0.0, 0.0, 0.0));
        let bv = b.vertex(Point3r::new(1.0, 0.0, 0.0));
        let c = b.vertex(Point3r::new(0.0, 1.0, 0.0));
        let d = b.vertex(Point3r::new(0.0, 0.0, 1.0));
        b.triangle(a, bv, c);
        b.triangle(a, c, d);
        b.triangle(a, d, bv);
        b.triangle(bv, d, c);
        b.build()
    }

    #[test]
    fn write_openfoam_creates_all_five_files() {
        let mesh = tet_mesh();
        let dir = std::env::temp_dir().join("cfd_mesh_of_test_five");
        write_openfoam_polymesh(&mesh, &dir, &[]).expect("write should succeed");

        for name in ["points", "faces", "owner", "neighbour", "boundary"] {
            assert!(
                dir.join(name).exists(),
                "expected file {name} to exist in output directory"
            );
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn points_file_has_correct_vertex_count() {
        let mesh = tet_mesh();
        let dir = std::env::temp_dir().join("cfd_mesh_of_test_vcount");
        write_openfoam_polymesh(&mesh, &dir, &[]).unwrap();
        let content = std::fs::read_to_string(dir.join("points")).unwrap();
        let count_line = content
            .lines()
            .find(|l| l.trim().parse::<usize>().is_ok())
            .unwrap_or("");
        assert_eq!(
            count_line.trim(),
            "4",
            "expected vertex count 4 in points file"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn faces_file_has_correct_face_count() {
        let mesh = tet_mesh();
        let dir = std::env::temp_dir().join("cfd_mesh_of_test_fcount");
        write_openfoam_polymesh(&mesh, &dir, &[]).unwrap();
        let content = std::fs::read_to_string(dir.join("faces")).unwrap();
        let count_line = content
            .lines()
            .filter(|l| l.trim().parse::<usize>().is_ok())
            .next()
            .unwrap_or("");
        assert_eq!(count_line.trim(), "4", "tetrahedron has 4 faces");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn boundary_file_default_patch_when_no_regions() {
        let mesh = tet_mesh();
        let dir = std::env::temp_dir().join("cfd_mesh_of_test_boundary");
        write_openfoam_polymesh(&mesh, &dir, &[]).unwrap();
        let content = std::fs::read_to_string(dir.join("boundary")).unwrap();
        assert!(
            content.contains("defaultFaces"),
            "boundary should contain 'defaultFaces' when no regions specified"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn boundary_file_named_patches_respected() {
        use crate::domain::mesh::IndexedMesh;

        // Build a mesh with two regions
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3r::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3r::new(0.0, 1.0, 0.0));
        let v3 = mesh.add_vertex_pos(Point3r::new(0.0, 0.0, 1.0));
        let r0 = RegionId::from_usize(0);
        let r1 = RegionId::from_usize(1);
        mesh.add_face_with_region(v0, v1, v2, r0);
        mesh.add_face_with_region(v0, v2, v3, r1);

        let dir = std::env::temp_dir().join("cfd_mesh_of_test_named");
        write_openfoam_polymesh(
            &mesh,
            &dir,
            &[
                (r0, "inlet", PatchType::Inlet),
                (r1, "outlet", PatchType::Outlet),
            ],
        )
        .unwrap();

        let content = std::fs::read_to_string(dir.join("boundary")).unwrap();
        assert!(content.contains("inlet"), "should contain 'inlet'");
        assert!(content.contains("outlet"), "should contain 'outlet'");
        assert!(
            content.contains("physicalType    inlet"),
            "inlet physicalType"
        );
        assert!(
            content.contains("physicalType    outlet"),
            "outlet physicalType"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn empty_mesh_returns_error() {
        let mesh = IndexedMesh::new();
        let dir = std::env::temp_dir().join("cfd_mesh_of_test_empty");
        let result = write_openfoam_polymesh(&mesh, &dir, &[]);
        assert!(result.is_err(), "empty mesh should produce an error");
    }
}
