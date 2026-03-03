"""Rewrite stitch_boundary_seams v4: merge-only, no face cleanup."""
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
    b"// boundary edges by merging their endpoint vertex IDs in the face set.\r\n"
    b"//\r\n"
    b"// IMPORTANT: this function does NOT remove degenerate or duplicate faces.\r\n"
    b"// That cleanup is deferred to `patch_small_boundary_holes` (Phase 6) which\r\n"
    b"// handles degenerate removal, duplicate removal, and non-manifold repair\r\n"
    b"// in the correct order to avoid topology cascade issues.\r\n"
    b"fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {\r\n"
    b"    // Build boundary edges: directed half-edges with no reverse partner.\r\n"
    b"    let mut he_count: HashMap<(VertexId, VertexId), u32> = HashMap::new();\r\n"
    b"    for face in faces.iter() {\r\n"
    b"        let v = face.vertices;\r\n"
    b"        for i in 0..3 {\r\n"
    b"            let j = (i + 1) % 3;\r\n"
    b"            *he_count.entry((v[i], v[j])).or_insert(0) += 1;\r\n"
    b"        }\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    let boundary_edges: Vec<(VertexId, VertexId)> = he_count\r\n"
    b"        .iter()\r\n"
    b"        .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he_count.contains_key(&(vj, vi)))\r\n"
    b"        .map(|(&e, _)| e)\r\n"
    b"        .collect();\r\n"
    b"\r\n"
    b"    if boundary_edges.is_empty() {\r\n"
    b"        return;\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    // Compute edge lengths.\r\n"
    b"    let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges\r\n"
    b"        .iter()\r\n"
    b"        .map(|&(vi, vj)| {\r\n"
    b"            let d = (pool.position(vj) - pool.position(vi)).norm_squared();\r\n"
    b"            (d, vi, vj)\r\n"
    b"        })\r\n"
    b"        .collect();\r\n"
    b"    edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));\r\n"
    b"\r\n"
    b"    let min_len_sq = edge_info.first().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"    let max_len_sq = edge_info.last().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"\r\n"
    b"    // Need a spread in edge lengths to identify short vs long boundary edges.\r\n"
    b"    if min_len_sq <= 0.0 || max_len_sq <= 0.0 || max_len_sq < 1.5 * min_len_sq {\r\n"
    b"        return;\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    // Threshold: geometric mean of min and max boundary edge lengths (squared).\r\n"
    b"    // For a bimodal distribution (short ~0.004, long ~0.05), this gives ~0.014.\r\n"
    b"    let threshold_sq = (min_len_sq * max_len_sq).sqrt();\r\n"
    b"\r\n"
    b"    // Build merge map: collapse every boundary edge shorter than threshold.\r\n"
    b"    let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"    for &(len_sq, vi, vj) in &edge_info {\r\n"
    b"        if len_sq >= threshold_sq {\r\n"
    b"            break; // sorted\r\n"
    b"        }\r\n"
    b"        // Chase existing merges to find canonical roots.\r\n"
    b"        let mut root_i = vi;\r\n"
    b"        while let Some(&next) = merge_map.get(&root_i) {\r\n"
    b"            root_i = next;\r\n"
    b"        }\r\n"
    b"        let mut root_j = vj;\r\n"
    b"        while let Some(&next) = merge_map.get(&root_j) {\r\n"
    b"            root_j = next;\r\n"
    b"        }\r\n"
    b"        if root_i == root_j {\r\n"
    b"            continue;\r\n"
    b"        }\r\n"
    b"        // Merge higher ID into lower for determinism.\r\n"
    b"        let (keep, discard) = if root_i < root_j { (root_i, root_j) } else { (root_j, root_i) };\r\n"
    b"        merge_map.insert(discard, keep);\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    if merge_map.is_empty() {\r\n"
    b"        return;\r\n"
    b"    }\r\n"
    b"\r\n"
    b"    #[cfg(test)]\r\n"
    b"    eprintln!(\r\n"
    b'        "[stitch] {} bnd edges, {} short (< {:.6}), {} merges",\r\n'
    b"        boundary_edges.len(),\r\n"
    b"        edge_info.iter().filter(|e| e.0 < threshold_sq).count(),\r\n"
    b"        threshold_sq.sqrt(), merge_map.len(),\r\n"
    b"    );\r\n"
    b"\r\n"
    b"    // Apply merge map to all faces (vertex ID remapping only, no face removal).\r\n"
    b"    for face in faces.iter_mut() {\r\n"
    b"        for v in &mut face.vertices {\r\n"
    b"            let mut t = *v;\r\n"
    b"            while let Some(&next) = merge_map.get(&t) {\r\n"
    b"                t = next;\r\n"
    b"            }\r\n"
    b"            *v = t;\r\n"
    b"        }\r\n"
    b"    }\r\n"
    b"    // Degenerate/duplicate face removal is handled by patch_small_boundary_holes.\r\n"
    b"}\r\n"
    b"\r\n"
)

new_data = data[:comment_start] + new_fn + data[found_end:]

with open(path, 'wb') as f:
    f.write(new_data)

print(f"SUCCESS: wrote {len(new_data)} bytes (was {len(data)} bytes)")
