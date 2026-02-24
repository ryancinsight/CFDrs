"""Fix stale mesh API calls in cfd-3d fem files."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src\fem"
FILES = ["solver.rs", "projection_solver.rs"]

FID = "cfd_mesh::domain::core::index::FaceId"
VID = "cfd_mesh::domain::core::index::VertexId"

for fname in FILES:
    path = os.path.join(BASE, fname)
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()

    original = content

    # Fix mesh.vertex(expr).unwrap().position.coords
    # Captures any expression inside (idx, v_i, v_j, m_idx, fallback, etc.)
    content = re.sub(
        r'mesh\.vertex\(([^)]+)\)\.unwrap\(\)\.position\.coords',
        lambda m: f'mesh.vertices.position({VID}::from_usize({m.group(1)})).coords',
        content
    )

    # Fix mesh.vertices().iter() to mesh.vertices.iter()
    content = content.replace("mesh.vertices().iter()", "mesh.vertices.iter()")

    # Fix "for v in mesh.vertices()" to "for v in mesh.vertices.iter()"
    content = content.replace("for v in mesh.vertices()", "for v in mesh.vertices.iter()")

    # Fix mesh.face(f_idx) - match "if let Some(f) = mesh.face(f_idx) {"
    # Replace with bounds-checked access using FaceStore
    content = re.sub(
        r'if let Some\((\w+)\) = mesh\.face\((\w+)\) \{',
        lambda m: f'if {m.group(2)} < mesh.face_count() {{ let {m.group(1)} = mesh.faces.get({FID}::from_usize({m.group(2)})); {{',
        content
    )

    # Fix "mesh.vertices().iter()" in vertex_positions collection (line with iter().map)
    content = re.sub(
        r'problem\.mesh\.vertices\(\)\.iter\(\)',
        'problem.mesh.vertices.iter()',
        content
    )

    # Fix integer comparison issue: local_verts.len() >= 4 (T vs integer)
    # The struct has T = generic but compares with integer literal - this is actual Rust integer
    # The error "expected type parameter T, found integer" means the len() result compared to T
    # This is actually a different error from `>= 4` being fine - it's another occurrence
    # where the comparison was accidentally changed to expect T. Let's not change the logic.

    # Add FaceId import at top of files that use it if not already present
    if f'{FID}' in content and 'use cfd_mesh::domain::core::index::FaceId' not in content and 'use cfd_mesh::domain::core::index::{FaceId' not in content:
        # Add after the last "use cfd_mesh" line
        content = re.sub(
            r'(use cfd_mesh[^;]+;\n)',
            r'\1use cfd_mesh::domain::core::index::{FaceId, VertexId};\n',
            content,
            count=1
        )

    if content != original:
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Fixed: {fname}")
    else:
        print(f"No changes: {fname}")
