"""Rewrite stitch_boundary_seams v2: chain-aware cross-seam merge."""
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

# Find the end of the function
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

# Include trailing \r\n
if data[found_end:found_end+2] == b'\r\n':
    found_end += 2
if data[found_end:found_end+2] == b'\r\n':
    found_end += 2

print(f"Old function: bytes {comment_start}..{found_end}")

new_fn = (
    b"// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.\r\n"
    b"//\r\n"
    b"// When two surfaces meet at a shallow angle, the CDT co-refinement may produce\r\n"
    b"// boundary edges that zigzag between the two surfaces: short cross-seam edges\r\n"
    b"// alternate with longer along-seam edges.  This function identifies the two\r\n"
    b"// seam chains via long-edge-only connectivity, then merges vertex pairs that\r\n"
    b"// are connected by short cross-seam edges across different chains.\r\n"
    b"//\r\n"
    b"// Algorithm:\r\n"
    b"//   1. Build half-edge map, identify boundary edges (valence 1).\r\n"
    b"//   2. Compute bimodal threshold (geometric mean of min/max edge lengths).\r\n"
    b"//   3. Build connectivity graph from long edges only -> connected components = chains.\r\n"
    b"//   4. For each short cross-seam edge connecting different chains, merge endpoints.\r\n"
    b"//   5. Apply merge map, remove degenerate/duplicate faces.\r\n"
    b"//   6. Repeat until no more cross-seam short edges or no boundary edges remain.\r\n"
    b"fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {\r\n"
    b"    const MAX_STITCH_ITERS: usize = 6;\r\n"
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
    b"        // Step 2: compute edge lengths and bimodal threshold.\r\n"
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
    b"        // No bimodal split if all edges are roughly the same length.\r\n"
    b"        if min_len_sq <= 0.0 || max_len_sq <= 0.0 || max_len_sq < 4.0 * min_len_sq {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Threshold = geometric mean of min and max edge lengths (squared).\r\n"
    b"        let threshold_sq = (min_len_sq * max_len_sq).sqrt();\r\n"
    b"\r\n"
    b"        // Classify edges.\r\n"
    b"        let short_edges: Vec<(VertexId, VertexId)> = edge_info\r\n"
    b"            .iter()\r\n"
    b"            .filter(|e| e.0 < threshold_sq)\r\n"
    b"            .map(|e| (e.1, e.2))\r\n"
    b"            .collect();\r\n"
    b"        let long_edges: Vec<(VertexId, VertexId)> = edge_info\r\n"
    b"            .iter()\r\n"
    b"            .filter(|e| e.0 >= threshold_sq)\r\n"
    b"            .map(|e| (e.1, e.2))\r\n"
    b"            .collect();\r\n"
    b"\r\n"
    b"        if short_edges.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 3: build connected components from LONG edges only.\r\n"
    b"        // Each component is a chain on one surface.\r\n"
    b"        let mut long_adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();\r\n"
    b"        for &(vi, vj) in &long_edges {\r\n"
    b"            long_adj.entry(vi).or_default().push(vj);\r\n"
    b"            long_adj.entry(vj).or_default().push(vi);\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Also add isolated vertices from short edges that have no long edge.\r\n"
    b"        let mut all_bnd_verts: Vec<VertexId> = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .flat_map(|&(a, b)| [a, b])\r\n"
    b"            .collect();\r\n"
    b"        all_bnd_verts.sort();\r\n"
    b"        all_bnd_verts.dedup();\r\n"
    b"\r\n"
    b"        let mut chain_id: HashMap<VertexId, usize> = HashMap::new();\r\n"
    b"        let mut current_chain: usize = 0;\r\n"
    b"        for &v in &all_bnd_verts {\r\n"
    b"            if chain_id.contains_key(&v) {\r\n"
    b"                continue;\r\n"
    b"            }\r\n"
    b"            // BFS using long edges only.\r\n"
    b"            let mut queue = vec![v];\r\n"
    b"            while let Some(u) = queue.pop() {\r\n"
    b"                if chain_id.contains_key(&u) {\r\n"
    b"                    continue;\r\n"
    b"                }\r\n"
    b"                chain_id.insert(u, current_chain);\r\n"
    b"                if let Some(nbrs) = long_adj.get(&u) {\r\n"
    b"                    for &n in nbrs {\r\n"
    b"                        if !chain_id.contains_key(&n) {\r\n"
    b"                            queue.push(n);\r\n"
    b"                        }\r\n"
    b"                    }\r\n"
    b"                }\r\n"
    b"            }\r\n"
    b"            current_chain += 1;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 4: for each short edge connecting DIFFERENT chains, merge endpoints.\r\n"
    b"        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"        for &(vi, vj) in &short_edges {\r\n"
    b"            let ci = chain_id.get(&vi).copied().unwrap_or(usize::MAX);\r\n"
    b"            let cj = chain_id.get(&vj).copied().unwrap_or(usize::MAX);\r\n"
    b"            if ci == cj {\r\n"
    b"                continue; // same chain, don't collapse\r\n"
    b"            }\r\n"
    b"            // Chase existing merges.\r\n"
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
    b"        // Also merge short same-chain edges only if they are VERY short\r\n"
    b"        // (< 1/16 of the geometric mean threshold) -- these are residual\r\n"
    b"        // near-duplicate vertices from arithmetic drift.\r\n"
    b"        let micro_tol_sq = threshold_sq * 0.0625;\r\n"
    b"        for &(len_sq, vi, vj) in &edge_info {\r\n"
    b"            if len_sq >= micro_tol_sq {\r\n"
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
    b"            if root_i != root_j {\r\n"
    b"                let (keep, discard) = if root_i < root_j { (root_i, root_j) } else { (root_j, root_i) };\r\n"
    b"                merge_map.insert(discard, keep);\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        if merge_map.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        #[cfg(test)]\r\n"
    b"        {\r\n"
    b"            let cross_count = short_edges.iter().filter(|&&(vi, vj)| {\r\n"
    b"                chain_id.get(&vi) != chain_id.get(&vj)\r\n"
    b"            }).count();\r\n"
    b"            eprintln!(\r\n"
    b'                "[stitch-iter {}] {} bnd edges, {} short, {} cross-chain, {} chains, {} merges, threshold={:.6}",\r\n'
    b"                _iter, boundary_edges.len(), short_edges.len(), cross_count,\r\n"
    b"                current_chain, merge_map.len(), threshold_sq.sqrt(),\r\n"
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
    b"        // Remove degenerate faces.\r\n"
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
