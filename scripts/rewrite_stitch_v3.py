"""Rewrite stitch_boundary_seams v3: aggressive short-edge collapse + wide merge."""
import sys

path = 'C:/Users/RyanClanton/gcli/software/millifluidic_design/CFDrs/crates/cfd-mesh/src/application/csg/arrangement/mod.rs'
with open(path, 'rb') as f:
    data = f.read()

fn_start_marker = b'fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {'
fn_start = data.find(fn_start_marker)
if fn_start < 0:
    print("ERROR: could not find stitch_boundary_seams function")
    sys.exit(1)

comment_marker = b'// Phase 5.5 helper: stitch boundary seams'
comment_start = data.rfind(comment_marker, 0, fn_start)
if comment_start < 0:
    print("ERROR: could not find comment")
    sys.exit(1)

pos = fn_start
brace_depth = 0
found_end = -1
while pos < len(data):
    b = data[pos]
    if b == ord('{'):
        brace_depth += 1
    elif b == ord('}'):
        brace_depth -= 1
        if brace_depth == 0:
            found_end = pos + 1
            break
    pos += 1

if found_end < 0:
    print("ERROR: could not find end of function")
    sys.exit(1)

if data[found_end:found_end+2] == b'\r\n':
    found_end += 2
if data[found_end:found_end+2] == b'\r\n':
    found_end += 2

new_fn = (
    b"// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.\r\n"
    b"//\r\n"
    b"// When two surfaces meet at a shallow angle, the CDT co-refinement may produce\r\n"
    b"// boundary edges that zigzag between the two surfaces: short cross-seam edges\r\n"
    b"// alternate with longer along-seam edges.  This function collapses the short\r\n"
    b"// boundary edges to close the seam, then uses a wider-tolerance all-pairs\r\n"
    b"// merge to clean up any residual near-coincident boundary vertices.\r\n"
    b"fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {\r\n"
    b"    // Helper: build boundary edges from the current face set.\r\n"
    b"    let build_bnd = |faces: &[FaceData]| -> Vec<(VertexId, VertexId)> {\r\n"
    b"        let mut he: HashMap<(VertexId, VertexId), u32> = HashMap::new();\r\n"
    b"        for face in faces {\r\n"
    b"            let v = face.vertices;\r\n"
    b"            for i in 0..3 {\r\n"
    b"                let j = (i + 1) % 3;\r\n"
    b"                *he.entry((v[i], v[j])).or_insert(0) += 1;\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"        he.iter()\r\n"
    b"            .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he.contains_key(&(vj, vi)))\r\n"
    b"            .map(|(&e, _)| e)\r\n"
    b"            .collect()\r\n"
    b"    };\r\n"
    b"\r\n"
    b"    // Helper: apply merge map with chain chasing, remove degenerate + duplicate faces.\r\n"
    b"    let apply_merge = |faces: &mut Vec<FaceData>, merge: &HashMap<VertexId, VertexId>| {\r\n"
    b"        if merge.is_empty() {\r\n"
    b"            return;\r\n"
    b"        }\r\n"
    b"        for face in faces.iter_mut() {\r\n"
    b"            for v in &mut face.vertices {\r\n"
    b"                let mut t = *v;\r\n"
    b"                while let Some(&next) = merge.get(&t) {\r\n"
    b"                    t = next;\r\n"
    b"                }\r\n"
    b"                *v = t;\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"        faces.retain(|f| {\r\n"
    b"            f.vertices[0] != f.vertices[1]\r\n"
    b"                && f.vertices[1] != f.vertices[2]\r\n"
    b"                && f.vertices[2] != f.vertices[0]\r\n"
    b"        });\r\n"
    b"        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());\r\n"
    b"        faces.retain(|f| {\r\n"
    b"            let mut key = f.vertices;\r\n"
    b"            key.sort();\r\n"
    b"            seen.insert(key)\r\n"
    b"        });\r\n"
    b"    };\r\n"
    b"\r\n"
    b"    // === Pass 1: iterative short-edge collapse ===\r\n"
    b"    for _iter in 0..6_usize {\r\n"
    b"        let boundary_edges = build_bnd(faces);\r\n"
    b"        if boundary_edges.is_empty() {\r\n"
    b"            return;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Compute edge lengths.\r\n"
    b"        let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .map(|&(vi, vj)| {\r\n"
    b"                let d = (pool.position(vj) - pool.position(vi)).norm_squared();\r\n"
    b"                (d, vi, vj)\r\n"
    b"            })\r\n"
    b"            .collect();\r\n"
    b"        edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));\r\n"
    b"\r\n"
    b"        let min_len_sq = edge_info.first().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"        let max_len_sq = edge_info.last().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"\r\n"
    b"        // Need at least 2x spread to identify a bimodal distribution.\r\n"
    b"        if min_len_sq <= 0.0 || max_len_sq < 2.0 * min_len_sq {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Threshold: geometric mean of min and max.\r\n"
    b"        let threshold_sq = (min_len_sq * max_len_sq).sqrt();\r\n"
    b"\r\n"
    b"        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"        for &(len_sq, vi, vj) in &edge_info {\r\n"
    b"            if len_sq >= threshold_sq {\r\n"
    b"                break;\r\n"
    b"            }\r\n"
    b"            let mut root_i = vi;\r\n"
    b"            while let Some(&next) = merge_map.get(&root_i) {\r\n"
    b"                root_i = next;\r\n"
    b"            }\r\n"
    b"            let mut root_j = vj;\r\n"
    b"            while let Some(&next) = merge_map.get(&root_j) {\r\n"
    b"                root_j = next;\r\n"
    b"            }\r\n"
    b"            if root_i == root_j {\r\n"
    b"                continue;\r\n"
    b"            }\r\n"
    b"            let (keep, discard) = if root_i < root_j { (root_i, root_j) } else { (root_j, root_i) };\r\n"
    b"            merge_map.insert(discard, keep);\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        if merge_map.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        #[cfg(test)]\r\n"
    b"        eprintln!(\r\n"
    b'            "[stitch-p1 {}] {} bnd, {} short (< {:.6}), {} merges",\r\n'
    b"            _iter, boundary_edges.len(),\r\n"
    b"            edge_info.iter().filter(|e| e.0 < threshold_sq).count(),\r\n"
    b"            threshold_sq.sqrt(), merge_map.len(),\r\n"
    b"        );\r\n"
    b"\r\n"
    b"        apply_merge(faces, &merge_map);\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    // === Pass 2: wide-tolerance nearest-boundary-vertex pair merge ===\r\n"
    b"    // After short-edge collapse, residual boundary vertices may be slightly\r\n"
    b"    // further apart than the collapse threshold.  Merge any boundary vertex\r\n"
    b"    // pair within a generous adaptive tolerance.\r\n"
    b"    for _iter in 0..4_usize {\r\n"
    b"        let boundary_edges = build_bnd(faces);\r\n"
    b"        if boundary_edges.is_empty() {\r\n"
    b"            return;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        let mut bnd_verts: Vec<VertexId> = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .flat_map(|&(a, b)| [a, b])\r\n"
    b"            .collect();\r\n"
    b"        bnd_verts.sort();\r\n"
    b"        bnd_verts.dedup();\r\n"
    b"\r\n"
    b"        // Adaptive tolerance: 50% of the average boundary edge length.\r\n"
    b"        let avg_len_sq: Real = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .map(|&(vi, vj)| (pool.position(vj) - pool.position(vi)).norm_squared())\r\n"
    b"            .sum::<Real>()\r\n"
    b"            / boundary_edges.len().max(1) as Real;\r\n"
    b"        let wide_tol_sq = avg_len_sq * 0.25; // (0.5 * avg_len)^2\r\n"
    b"\r\n"
    b"        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"        for i in 0..bnd_verts.len() {\r\n"
    b"            let vi = bnd_verts[i];\r\n"
    b"            if merge_map.contains_key(&vi) {\r\n"
    b"                continue;\r\n"
    b"            }\r\n"
    b"            let pi = pool.position(vi);\r\n"
    b"            let mut best_d = wide_tol_sq;\r\n"
    b"            let mut best_j: Option<VertexId> = None;\r\n"
    b"            for j in (i + 1)..bnd_verts.len() {\r\n"
    b"                let vj = bnd_verts[j];\r\n"
    b"                if merge_map.contains_key(&vj) {\r\n"
    b"                    continue;\r\n"
    b"                }\r\n"
    b"                let d = (pool.position(vj) - pi).norm_squared();\r\n"
    b"                if d < best_d {\r\n"
    b"                    best_d = d;\r\n"
    b"                    best_j = Some(vj);\r\n"
    b"                }\r\n"
    b"            }\r\n"
    b"            if let Some(vj) = best_j {\r\n"
    b"                merge_map.insert(vj, vi); // vi < vj since sorted\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        if merge_map.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        #[cfg(test)]\r\n"
    b"        eprintln!(\r\n"
    b'            "[stitch-p2 {}] {} bnd edges, {} bnd verts, tol={:.6}, {} merges",\r\n'
    b"            _iter, boundary_edges.len(), bnd_verts.len(),\r\n"
    b"            wide_tol_sq.sqrt(), merge_map.len(),\r\n"
    b"        );\r\n"
    b"\r\n"
    b"        apply_merge(faces, &merge_map);\r\n"
    b"    }\r\n"
    b"}\r\n"
    b"\r\n"
)

new_data = data[:comment_start] + new_fn + data[found_end:]

with open(path, 'wb') as f:
    f.write(new_data)

print(f"SUCCESS: wrote {len(new_data)} bytes (was {len(data)} bytes)")
