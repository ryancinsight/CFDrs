"""Insert Phase 5.5 stitch_boundary_seams into mod.rs."""
import sys

path = 'C:/Users/RyanClanton/gcli/software/millifluidic_design/CFDrs/crates/cfd-mesh/src/application/csg/arrangement/mod.rs'
with open(path, 'rb') as f:
    data = f.read()

# Find insertion point: right after "result_faces.extend(kept_faces);\n\r\n"
target = b'result_faces.extend(kept_faces);'
idx = data.find(target)
if idx < 0:
    print("ERROR: could not find target line")
    sys.exit(1)
eol = data.index(b'\n', idx)
insert_after = eol + 1  # after \n
# skip \r\n blank line
if data[insert_after:insert_after+2] == b'\r\n':
    insert_after += 2

# The Phase 5.5 call to insert
phase55_call = (
    b"\r\n"
    b"    // Phase 5.5: stitch boundary seams from unresolved intersection gaps.\r\n"
    b"    // At shallow angles, some face pairs barely intersect and CDT co-refinement\r\n"
    b"    // produces separate seam chains ~10-20um apart.  This pass finds and merges\r\n"
    b"    // corresponding boundary vertex pairs to close those seams.\r\n"
    b"    stitch_boundary_seams(&mut result_faces, pool);\r\n"
    b"\r\n"
)

# Find where to insert the function definition (before patch_small_boundary_holes).
fn_target = b'fn patch_small_boundary_holes(faces: &mut Vec<FaceData>, pool: &VertexPool) {'
fn_idx = data.find(fn_target)
if fn_idx < 0:
    print("ERROR: could not find patch_small_boundary_holes function")
    sys.exit(1)

# Go back to the "Phase 6 helper" comment line
phase6_helper = data.rfind(b'Phase 6 helper', 0, fn_idx)
if phase6_helper >= 0:
    line_start = data.rfind(b'\n', 0, phase6_helper)
    fn_insert = line_start + 1 if line_start >= 0 else phase6_helper
else:
    # Fallback: go back from fn_idx to find start of line
    line_start = data.rfind(b'\n', 0, fn_idx)
    fn_insert = line_start + 1 if line_start >= 0 else fn_idx

print(f"Phase 5.5 call insert at byte {insert_after}")
print(f"Function definition insert at byte {fn_insert}")

# Build the function definition
fn_def = (
    b"// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.\r\n"
    b"//\r\n"
    b"// When two surfaces meet at a shallow angle, the CDT co-refinement may produce\r\n"
    b"// separate boundary chains on opposite sides of the intersection seam, with\r\n"
    b"// vertices 10-20um apart.  This function detects such paired chains and merges\r\n"
    b"// corresponding boundary vertices to produce manifold edges.\r\n"
    b"//\r\n"
    b"// Algorithm:\r\n"
    b"//   1. Build half-edge map, identify boundary edges (valence 1).\r\n"
    b"//   2. Collect boundary vertices, build spatial hash.\r\n"
    b"//   3. Compute adaptive tolerance from min boundary edge length.\r\n"
    b"//   4. For each boundary vertex, find nearest boundary vertex within tolerance\r\n"
    b"//      that is NOT on the same chain (connected via boundary edges).\r\n"
    b"//   5. Merge matched pairs, remove degenerate/duplicate faces.\r\n"
    b"//   6. Repeat up to 3 iterations until boundary count stabilizes.\r\n"
    b"fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {\r\n"
    b"    const MAX_STITCH_ITERS: usize = 3;\r\n"
    b"\r\n"
    b"    for iter in 0..MAX_STITCH_ITERS {\r\n"
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
    b"        // Step 2: collect unique boundary vertices.\r\n"
    b"        let mut bnd_verts: Vec<VertexId> = boundary_edges\r\n"
    b"            .iter()\r\n"
    b"            .flat_map(|&(a, b)| [a, b])\r\n"
    b"            .collect();\r\n"
    b"        bnd_verts.sort();\r\n"
    b"        bnd_verts.dedup();\r\n"
    b"\r\n"
    b"        if bnd_verts.is_empty() {\r\n"
    b"            break;\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 3: compute adaptive tolerance from min boundary edge length.\r\n"
    b"        // Use 25% of the minimum boundary edge length, capped at 50um (2.5e-9 sq).\r\n"
    b"        let mut min_edge_len_sq: Real = Real::MAX;\r\n"
    b"        for &(vi, vj) in &boundary_edges {\r\n"
    b"            let pi = pool.position(vi);\r\n"
    b"            let pj = pool.position(vj);\r\n"
    b"            let d = (pj - pi).norm_squared();\r\n"
    b"            if d > 1e-20 && d < min_edge_len_sq {\r\n"
    b"                min_edge_len_sq = d;\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"        // Adaptive tolerance: 25% of min edge length, squared.\r\n"
    b"        // For typical 14um gaps with ~0.5mm edges, this gives ~125um tolerance.\r\n"
    b"        // Cap at 50um^2 = 2.5e-9 m^2 to avoid over-merging.\r\n"
    b"        let adaptive_tol_sq = (0.0625 * min_edge_len_sq).min(2.5e-9);\r\n"
    b"\r\n"
    b"        // Step 4: build adjacency for chain membership (BFS connected components).\r\n"
    b"        let mut bnd_adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();\r\n"
    b"        for &(vi, vj) in &boundary_edges {\r\n"
    b"            bnd_adj.entry(vi).or_default().push(vj);\r\n"
    b"            bnd_adj.entry(vj).or_default().push(vi);\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Assign chain IDs via BFS.\r\n"
    b"        let mut chain_id: HashMap<VertexId, usize> = HashMap::new();\r\n"
    b"        let mut current_chain: usize = 0;\r\n"
    b"        for &v in &bnd_verts {\r\n"
    b"            if chain_id.contains_key(&v) {\r\n"
    b"                continue;\r\n"
    b"            }\r\n"
    b"            let mut queue = vec![v];\r\n"
    b"            while let Some(u) = queue.pop() {\r\n"
    b"                if chain_id.contains_key(&u) {\r\n"
    b"                    continue;\r\n"
    b"                }\r\n"
    b"                chain_id.insert(u, current_chain);\r\n"
    b"                if let Some(nbrs) = bnd_adj.get(&u) {\r\n"
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
    b"        // Step 5: for each boundary vertex, find nearest vertex on a DIFFERENT chain.\r\n"
    b"        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();\r\n"
    b"        for i in 0..bnd_verts.len() {\r\n"
    b"            let vi = bnd_verts[i];\r\n"
    b"            if merge_map.contains_key(&vi) {\r\n"
    b"                continue;\r\n"
    b"            }\r\n"
    b"            let ci = chain_id[&vi];\r\n"
    b"            let pi = pool.position(vi);\r\n"
    b"            let mut best_dist_sq: Real = adaptive_tol_sq;\r\n"
    b"            let mut best_vj: Option<VertexId> = None;\r\n"
    b"\r\n"
    b"            for j in 0..bnd_verts.len() {\r\n"
    b"                if i == j {\r\n"
    b"                    continue;\r\n"
    b"                }\r\n"
    b"                let vj = bnd_verts[j];\r\n"
    b"                if merge_map.contains_key(&vj) {\r\n"
    b"                    continue;\r\n"
    b"                }\r\n"
    b"                let cj = chain_id[&vj];\r\n"
    b"                if ci == cj {\r\n"
    b"                    continue; // same chain -- skip\r\n"
    b"                }\r\n"
    b"                let pj = pool.position(vj);\r\n"
    b"                let d = (pj - pi).norm_squared();\r\n"
    b"                if d < best_dist_sq {\r\n"
    b"                    best_dist_sq = d;\r\n"
    b"                    best_vj = Some(vj);\r\n"
    b"                }\r\n"
    b"            }\r\n"
    b"\r\n"
    b"            if let Some(vj) = best_vj {\r\n"
    b"                // Merge higher ID into lower ID for determinism.\r\n"
    b"                let (keep, discard) = if vi < vj { (vi, vj) } else { (vj, vi) };\r\n"
    b"                merge_map.insert(discard, keep);\r\n"
    b"            }\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        if merge_map.is_empty() {\r\n"
    b"            break; // no progress -- stop iterating\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        #[cfg(test)]\r\n"
    b"        {\r\n"
    b"            eprintln!(\r\n"
    b'                "[stitch-iter {}] {} boundary edges, {} merge pairs, tol_sq={:.2e}",\r\n'
    b"                iter, boundary_edges.len(), merge_map.len(), adaptive_tol_sq\r\n"
    b"            );\r\n"
    b"        }\r\n"
    b"\r\n"
    b"        // Step 6: apply merge map to all faces, remove degenerate/duplicate.\r\n"
    b"        for face in faces.iter_mut() {\r\n"
    b"            for v in &mut face.vertices {\r\n"
    b"                let mut t = *v;\r\n"
    b"                while let Some(&next) = merge_map.get(&t) {\r\n"
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

# Perform the two insertions (higher offset first to avoid shifting).
assert fn_insert > insert_after, "Expected fn_insert > insert_after"
new_data = data[:fn_insert] + fn_def + data[fn_insert:]
new_data = new_data[:insert_after] + phase55_call + new_data[insert_after:]

with open(path, 'wb') as f:
    f.write(new_data)

print(f"SUCCESS: wrote {len(new_data)} bytes (was {len(data)})")
print(f"Inserted {len(phase55_call)} bytes for Phase 5.5 call")
print(f"Inserted {len(fn_def)} bytes for function definition")
