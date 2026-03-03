"""Rewrite stitch_boundary_seams to use short-edge-collapse instead of chain matching."""
import sys

path = 'C:/Users/RyanClanton/gcli/software/millifluidic_design/CFDrs/crates/cfd-mesh/src/application/csg/arrangement/mod.rs'
with open(path, 'rb') as f:
    data = f.read()

# Find the old function
fn_start_marker = b'fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {'
fn_start = data.find(fn_start_marker)
if fn_start < 0:
    print("ERROR: could not find stitch_boundary_seams function")
    sys.exit(1)

# Find the comment block above the function
comment_marker = b'// Phase 5.5 helper: stitch boundary seams'
comment_start = data.rfind(comment_marker, 0, fn_start)
if comment_start < 0:
    print("ERROR: could not find Phase 5.5 comment")
    sys.exit(1)

# Find the end of the function - look for the closing brace at the same indent level
# The function ends with \r\n}\r\n at the top level (no indentation before })
# Search for \n}\r\n after fn_start
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
            found_end = pos + 1  # include the }
            break
    pos += 1

if found_end < 0:
    print("ERROR: could not find end of stitch_boundary_seams function")
    sys.exit(1)

# Include trailing \r\n after the closing }
if data[found_end:found_end+2] == b'\r\n':
    found_end += 2
if data[found_end:found_end+2] == b'\r\n':
    found_end += 2

print(f"Old function: bytes {comment_start}..{found_end} ({found_end - comment_start} bytes)")
print(f"Context before: {repr(data[comment_start-20:comment_start])}")
print(f"Context after:  {repr(data[found_end:found_end+60])}")

new_fn = (
    b"// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.\r\n"
    b"//\r\n"
    b"// When two surfaces meet at a shallow angle, the CDT co-refinement may produce\r\n"
    b"// boundary edges that zigzag between the two surfaces: short cross-seam edges\r\n"
    b"// (~0.004-0.01 units) alternate with longer along-seam edges (~0.05 units).\r\n"
    b"// This function collapses the short boundary edges to close the seam.\r\n"
    b"//\r\n"
    b"// Algorithm:\r\n"
    b"//   1. Build half-edge map, identify boundary edges (valence 1).\r\n"
    b"//   2. Sort boundary edges by length.\r\n"
    b"//   3. Compute collapse threshold from the bimodal edge-length distribution:\r\n"
    b"//      geometric mean of min and max boundary edge lengths.\r\n"
    b"//   4. Collapse all boundary edges shorter than threshold.\r\n"
    b"//   5. Remove degenerate/duplicate faces.\r\n"
    b"//   6. Repeat up to 4 iterations until no short edges remain.\r\n"
    b"fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {\r\n"
    b"    const MAX_STITCH_ITERS: usize = 4;\r\n"
    b"\r\n"
    b"    for _iter in 0..MAX_STITCH_ITERS {\r\n"
    b"        // Step 1: build half-edge map and find boundary edges.\r\n"
    b"        let mut he_count: HashMap<(VertexId, VertexId), u32> = HashMap::new();\r\n"
    b"        for face in faces.iter() {\r\n"
    b"            let v = face.vertices;\r\n"
    b"            for i in 0..3 {\r\n"
    b"                let j = (i + 1) % 3;\r\n"
    b"                *he_count.entry((v[i], v[j])).or_insert(0) += 1;\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        let boundary_edges: Vec<(VertexId, VertexId)> = he_count\r\n"
    b"            .iter()\r\n"
    b"            .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he_count.contains_key(&(vj, vi)))\r\n"
    b"            .map(|(&e, _)| e)\r\n"
    b"            .collect();\r\n"
    b"\r\n"
    b"        if boundary_edges.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 2: compute edge lengths.\r\n"
    b"        let mut edge_lens_sq: Vec<(Real, VertexId, VertexId)> = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .map(|&(vi, vj)| {\r\n"
    b"                let d = (pool.position(vj) - pool.position(vi)).norm_squared();\r\n"
    b"                (d, vi, vj)\r\n"
    b"            })\r\n"
    b"            .collect();\r\n"
    b"        edge_lens_sq.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));\r\n"
    b"\r\n"
    b"        let min_len_sq = edge_lens_sq.first().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"        let max_len_sq = edge_lens_sq.last().map(|e| e.0).unwrap_or(0.0);\r\n"
    b"\r\n"
    b"        // If all edges are about the same length, no bimodal split to exploit.\r\n"
    b"        if min_len_sq <= 0.0 || max_len_sq <= 0.0 || max_len_sq < 4.0 * min_len_sq {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 3: collapse threshold = geometric mean of min and max edge lengths (squared).\r\n"
    b"        // For bimodal distribution (short ~0.004, long ~0.05), this gives ~0.014,\r\n"
    b"        // which cleanly separates the two populations.\r\n"
    b"        let threshold_sq = (min_len_sq * max_len_sq).sqrt();\r\n"
    b"\r\n"
    b"        // Step 4: collect vertex merge pairs from short edges.\r\n"
    b"        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"        for &(len_sq, vi, vj) in &edge_lens_sq {\r\n"
    b"            if len_sq >= threshold_sq {\r\n"
    b"                break; // sorted, so all remaining are longer\r\n"
    b"            }\r\n"
    b"            // Chase existing merges to find canonical vertex.\r\n"
    b"            let mut root_i = vi;\r\n"
    b"            while let Some(&next) = merge_map.get(&root_i) {\r\n"
    b"                root_i = next;\r\n"
    b"            }\r\n"
    b"            let mut root_j = vj;\r\n"
    b"            while let Some(&next) = merge_map.get(&root_j) {\r\n"
    b"                root_j = next;\r\n"
    b"            }\r\n"
    b"            if root_i == root_j {\r\n"
    b"                continue; // already merged\r\n"
    b"            }\r\n"
    b"            // Merge higher ID into lower ID for determinism.\r\n"
    b"            let (keep, discard) = if root_i < root_j { (root_i, root_j) } else { (root_j, root_i) };\r\n"
    b"            merge_map.insert(discard, keep);\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        if merge_map.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        #[cfg(test)]\r\n"
    b"        {\r\n"
    b"            eprintln!(\r\n"
    b'                "[stitch-iter {}] {} boundary edges, {} short (< {:.6}), {} merges",\r\n'
    b"                _iter,\r\n"
    b"                boundary_edges.len(),\r\n"
    b"                edge_lens_sq.iter().filter(|e| e.0 < threshold_sq).count(),\r\n"
    b"                threshold_sq.sqrt(),\r\n"
    b"                merge_map.len(),\r\n"
    b"            );\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 5: apply merge map to all faces.\r\n"
    b"        for face in faces.iter_mut() {\r\n"
    b"            for v in &mut face.vertices {\r\n"
    b"                let mut t = *v;\r\n"
    b"                while let Some(&next) = merge_map.get(&t) {\r\n"
    b"                    t = next;\r\n"
    b"                }\r\n"
    b"                *v = t;\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Remove degenerate faces (collapsed triangles).\r\n"
    b"        faces.retain(|f| {\r\n"
    b"            f.vertices[0] != f.vertices[1]\r\n"
    b"                && f.vertices[1] != f.vertices[2]\r\n"
    b"                && f.vertices[2] != f.vertices[0]\r\n"
    b"        });\r\n"
    b"\r\n"
    b"        // Remove duplicate faces.\r\n"
    b"        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());\r\n"
    b"        faces.retain(|f| {\r\n"
    b"            let mut key = f.vertices;\r\n"
    b"            key.sort();\r\n"
    b"            seen.insert(key)\r\n"
    b"        });\r\n"
    b"    }\r\n"
    b"}\r\n"
    b"\r\n"
)

new_data = data[:comment_start] + new_fn + data[found_end:]

with open(path, 'wb') as f:
    f.write(new_data)

print(f"SUCCESS: wrote {len(new_data)} bytes (was {len(data)} bytes)")
print(f"Replaced {found_end - comment_start} bytes with {len(new_fn)} bytes")
