import os

path = r'c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-mesh\src\domain\mesh.rs'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

mesh_dir = r'c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-mesh\src\domain\mesh'
os.makedirs(mesh_dir, exist_ok=True)

with open(os.path.join(mesh_dir, 'generic.rs'), 'w', encoding='utf-8') as f:
    f.writelines(lines[44:288])

with open(os.path.join(mesh_dir, 'indexed.rs'), 'w', encoding='utf-8') as f:
    f.write('use nalgebra::{Point3, Vector3};\n')
    f.writelines(lines[288:579])

with open(os.path.join(mesh_dir, 'halfedge.rs'), 'w', encoding='utf-8') as f:
    f.write('use nalgebra::Point3;\n')
    f.writelines(lines[579:])

mod_content = """//! # Mesh Types
//!
//! This module provides three mesh types:
//!
//! - **[`HalfEdgeMesh<'id>`]** — the new state-of-the-art mesh backed by a
//!   GhostCell-permissioned half-edge topology kernel with SlotMap generational
//!   keys. Use [`with_mesh`] as the entry point.
//! - **`Mesh<T>`** — legacy FEM/FVM mesh with typed vertices, faces, and cells.
//!   Kept for downstream compatibility with `cfd-3d` geometry builders.
//! - **[`IndexedMesh`]** — the watertight-first surface mesh combining
//!   `VertexPool` (spatial-hash dedup), `FaceStore`, `EdgeStore`, and
//!   `AttributeStore`. Also kept for backward compatibility.

pub mod generic;
pub mod halfedge;
pub mod indexed;

#[allow(deprecated)]
pub use generic::{Mesh, MeshStatistics};
pub use halfedge::{with_mesh, HalfEdgeMesh};
pub use indexed::{IndexedMesh, MeshBuilder};
"""

with open(os.path.join(mesh_dir, 'mod.rs'), 'w', encoding='utf-8') as f:
    f.write(mod_content)

os.remove(path)
