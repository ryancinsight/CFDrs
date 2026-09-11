//! AMG coarsening edge-case tests.
//!
//! These tests were originally parked empty "pending migration of the
//! domain-specific multigrid code to the leto-ops API surface". That
//! migration is complete: the coarsening suite is cfd-math-owned
//! (`cfd_math::multigrid`, scalar traits via leto), so the two edge
//! cases — undecided points after the first Ruge-Stüben pass, and
//! disconnected island components — are exercised for real
//! (CFDRS-GA-013).
//!
//! # Schemas under test
//!
//! - `CoarseningResult::fine_to_coarse_map[i] = Some(c)` requires
//!   `c < coarse_points.len()`, and C-points self-map
//!   (`coarse_points[c] == i`) — the mapping schema of the coarsening
//!   module's own correctness tests.
//! - Disconnected islands must never leak F→C associations across
//!   components: every mapped F-point's coarse target lies in the same
//!   connected component as the F-point.

use cfd_math::multigrid::{CoarseningResult, ruge_stueben_coarsening};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};

/// 5×5 two-dimensional five-point Laplacian (25 strongly connected
/// points) — the same schema as the coarsening module's unit helpers,
/// rebuilt through the public `SparseMatrixBuilder`.
fn grid_laplacian(n: usize) -> SparseMatrix<f64> {
    let size = n * n;
    let mut builder = SparseMatrixBuilder::new(size, size);
    for i in 0..n {
        for j in 0..n {
            let idx = i * n + j;
            if i > 0 {
                builder
                    .add_entry(idx, (i - 1) * n + j, -1.0)
                    .expect("invariant: grid neighbour indices are in range");
            }
            if j > 0 {
                builder
                    .add_entry(idx, i * n + (j - 1), -1.0)
                    .expect("invariant: grid neighbour indices are in range");
            }
            builder
                .add_entry(idx, idx, 4.0)
                .expect("invariant: grid diagonal index is in range");
            if j < n - 1 {
                builder
                    .add_entry(idx, i * n + (j + 1), -1.0)
                    .expect("invariant: grid neighbour indices are in range");
            }
            if i < n - 1 {
                builder
                    .add_entry(idx, (i + 1) * n + j, -1.0)
                    .expect("invariant: grid neighbour indices are in range");
            }
        }
    }
    builder
        .build()
        .expect("invariant: the assembled grid Laplacian is a valid CSR matrix")
}

/// 6-point block-diagonal matrix holding two disconnected 3-point
/// 1D path Laplacians — the "island" topology. No entry couples the
/// two blocks.
fn island_matrix() -> SparseMatrix<f64> {
    let points = 6;
    let block = 3;
    let mut builder = SparseMatrixBuilder::new(points, points);
    for island in 0..2 {
        let base = island * block;
        for k in 0..block {
            builder
                .add_entry(base + k, base + k, 2.0)
                .expect("invariant: island diagonal index is in range");
            if k > 0 {
                builder
                    .add_entry(base + k, base + k - 1, -1.0)
                    .expect("invariant: island neighbour index is in range");
            }
            if k < block - 1 {
                builder
                    .add_entry(base + k, base + k + 1, -1.0)
                    .expect("invariant: island neighbour index is in range");
            }
        }
    }
    builder
        .build()
        .expect("invariant: the assembled block-diagonal matrix is a valid CSR matrix")
}

/// After the first Ruge-Stüben pass the second pass classifies every
/// remaining undecided point as a C-point, so on termination the split
/// must be total: every point is either a C-point or an F-point with a
/// valid coarse mapping. A dangling undecided point (map `None` while
/// absent from `coarse_points`) would silently drop a DOF from the
/// coarse grid and corrupt the hierarchy.
#[test]
fn test_undecided_points_fully_classified() {
    let matrix = grid_laplacian(5);
    let result: CoarseningResult<f64> =
        ruge_stueben_coarsening(&matrix, 0.25).expect("invariant: coarsening inputs are valid");

    let n = matrix.nrows();
    assert_eq!(result.fine_to_coarse_map.len(), n);

    let mut classified = vec![false; n];
    for &c in &result.coarse_points {
        assert!(c < n, "coarse point index out of range");
        classified[c] = true;
    }

    // Mapping schema: a mapped point either self-maps (it is a C-point
    // and `coarse_points[c] == i`) or is an F-point pointing at a
    // distinct coarse representative. `None` is only legal for a
    // C-point — the second pass classifies every leftover undecided
    // point, so a `None` on a non-C-point would be a dropped DOF.
    for (i, map) in result.fine_to_coarse_map.iter().enumerate() {
        match map {
            Some(c) => {
                assert!(
                    *c < result.coarse_points.len(),
                    "F-point {i} maps out of the coarse set"
                );
                let representative = result.coarse_points[*c];
                if representative == i {
                    assert!(classified[i], "self-mapped point {i} must be a C-point");
                } else {
                    assert!(!classified[i], "mapped point {i} must be an F-point");
                }
            }
            None => {
                assert!(
                    classified[i],
                    "point {i} is undecided: neither C-point nor mapped F-point"
                );
            }
        }
    }

    // Every C-point self-maps at its own coarse index.
    for (idx, &cp) in result.coarse_points.iter().enumerate() {
        assert_eq!(
            result.fine_to_coarse_map[cp],
            Some(idx),
            "coarse point {cp} not correctly self-mapped"
        );
    }

    assert!(
        !result.coarse_points.is_empty(),
        "a connected grid must produce at least one C-point"
    );
}

/// Two disconnected islands must coarsen without cross-component
/// leakage: each island contributes at least one C-point, every F-point
/// maps to a coarse target inside its own island, and the union covers
/// all points. A cross-island association would interpolate from values
/// on a decoupled subdomain.
#[test]
fn test_island_points_stay_within_components() {
    let matrix = island_matrix();
    let result: CoarseningResult<f64> =
        ruge_stueben_coarsening(&matrix, 0.25).expect("invariant: coarsening inputs are valid");

    let n = matrix.nrows();
    let block = 3;
    let island_of = |p: usize| p / block;

    let mut island_has_c = [false; 2];
    for &c in &result.coarse_points {
        island_has_c[island_of(c)] = true;
    }
    assert!(
        island_has_c.iter().all(|&has| has),
        "every island must retain at least one C-point"
    );

    for (i, map) in result.fine_to_coarse_map.iter().enumerate() {
        if let Some(c) = map {
            let target = result.coarse_points[*c];
            assert_eq!(
                island_of(target),
                island_of(i),
                "F-point {i} in island {} maps to C-point {target} in island {}",
                island_of(i),
                island_of(target)
            );
        }
    }

    // Split totality: every point is either a C-point or a mapped
    // F-point — no DOF may be dropped from the coarse hierarchy.
    let mut classified = vec![false; n];
    for &c in &result.coarse_points {
        classified[c] = true;
    }
    for (i, map) in result.fine_to_coarse_map.iter().enumerate() {
        assert!(
            classified[i] || map.is_some(),
            "point {i} is neither C-point nor mapped F-point"
        );
    }
}
