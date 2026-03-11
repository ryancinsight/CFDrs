//! Shared fragment-analysis helpers for arrangement CSG.
//!
//! These utilities consolidate source-preserving component detection and
//! sliver rejection across the binary and generalized Boolean paths.

#[cfg(test)]
use super::classify::tri_normal;
use super::dsu::DisjointSet;
use crate::domain::core::constants::SLIVER_AREA_RATIO_SQ;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Vector3r};

/// Build connected-component roots for fragment adjacency induced by shared
/// edges, while preserving source boundaries.
pub(crate) fn component_roots_by_source<T, FVertices, FSource, S>(
    fragments: &[T],
    mut face_vertices: FVertices,
    mut source_key: FSource,
) -> Vec<usize>
where
    FVertices: FnMut(&T) -> [VertexId; 3],
    FSource: FnMut(&T) -> S,
    S: Copy + Eq,
{
    let mut edge_refs: Vec<(VertexId, VertexId, usize, S)> =
        Vec::with_capacity(fragments.len() * 3);
    for (fragment_index, fragment) in fragments.iter().enumerate() {
        let vertices = face_vertices(fragment);
        let source = source_key(fragment);
        for edge_index in 0..3 {
            let a = vertices[edge_index];
            let b = vertices[(edge_index + 1) % 3];
            let (mn, mx) = if a < b { (a, b) } else { (b, a) };
            edge_refs.push((mn, mx, fragment_index, source));
        }
    }
    edge_refs.sort_unstable_by_key(|&(a, b, _, _)| (a.raw(), b.raw()));

    let mut dsu = DisjointSet::new(fragments.len());
    let mut run_start = 0usize;
    while run_start < edge_refs.len() {
        let (edge_a, edge_b, _, first_source) = edge_refs[run_start];
        let mut run_end = run_start + 1;
        while run_end < edge_refs.len()
            && edge_refs[run_end].0 == edge_a
            && edge_refs[run_end].1 == edge_b
        {
            run_end += 1;
        }

        let mixed_sources = edge_refs[run_start..run_end]
            .iter()
            .any(|&(_, _, _, source)| source != first_source);
        if mixed_sources {
            run_start = run_end;
            continue;
        }

        if run_end - run_start >= 2 {
            let root = edge_refs[run_start].2;
            for &(_, _, index, _) in &edge_refs[(run_start + 1)..run_end] {
                dsu.union(root, index);
            }
        }
        run_start = run_end;
    }

    let mut roots = Vec::with_capacity(fragments.len());
    for index in 0..fragments.len() {
        roots.push(dsu.find(index));
    }
    roots
}

/// Return whether a triangle is a numerically degenerate sliver under the
/// canonical arrangement threshold.
pub(crate) fn is_degenerate_sliver_with_normal(tri: &[Point3r; 3], normal: &Vector3r) -> bool {
    let p0 = tri[0];
    let p1 = tri[1];
    let p2 = tri[2];

    let e01_sq = (p1 - p0).norm_squared();
    let e02_sq = (p2 - p0).norm_squared();
    let e12_sq = (p2 - p1).norm_squared();
    let area_sq = normal.norm_squared();
    let max_edge_sq = e01_sq.max(e02_sq).max(e12_sq);

    max_edge_sq > 1e-20 && area_sq < SLIVER_AREA_RATIO_SQ * max_edge_sq
}

/// Convenience wrapper that computes the normal before evaluating the sliver
/// criterion.
#[cfg(test)]
pub(crate) fn is_degenerate_sliver(tri: &[Point3r; 3]) -> bool {
    let normal = tri_normal(tri);
    is_degenerate_sliver_with_normal(tri, &normal)
}
