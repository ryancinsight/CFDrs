//! Tests for MPI functionality

use super::decomposition::{
    DecompositionStrategy, DomainDecomposition, NeighborDirection, NeighborInfo,
};
use super::distributed_solvers::{
    AdditiveSchwarzPreconditioner, BlockJacobiPreconditioner, DistributedGMRES,
    DistributedLinearOperator,
};
use super::*;
use nalgebra::RealField;

#[test]
fn test_global_extents_creation() {
    // Test 2D global extents
    let extents_2d = GlobalExtents::new_2d(100, 50, (0.0, 1.0, 0.0, 0.5));
    assert_eq!(extents_2d.nx_global, 100);
    assert_eq!(extents_2d.ny_global, 50);
    assert_eq!(extents_2d.nz_global, 1);
    assert_eq!(extents_2d.bounds, (0.0, 1.0, 0.0, 0.5, 0.0, 0.0));

    // Test 3D global extents
    let extents_3d = GlobalExtents::new_3d(50, 40, 30, (0.0, 2.0, 0.0, 1.5, 0.0, 3.0));
    assert_eq!(extents_3d.nx_global, 50);
    assert_eq!(extents_3d.ny_global, 40);
    assert_eq!(extents_3d.nz_global, 30);
    assert_eq!(extents_3d.bounds, (0.0, 2.0, 0.0, 1.5, 0.0, 3.0));
}

#[test]
fn test_local_subdomain_operations() {
    let subdomain = LocalSubdomain {
        rank: 1,
        nx_local: 25,
        ny_local: 20,
        nz_local: 1,
        i_start_global: 25,
        j_start_global: 0,
        k_start_global: 0,
        ghost_layers: 1,
    };

    // Test total dimensions with ghosts
    assert_eq!(subdomain.total_nx(), 27); // 25 + 2*1
    assert_eq!(subdomain.total_ny(), 22); // 20 + 2*1
    assert_eq!(subdomain.total_nz(), 3); // 1 + 2*1

    // Test owned cell checks
    assert!(subdomain.is_owned_cell(1, 1, 0)); // First owned cell
    assert!(subdomain.is_owned_cell(25, 19, 0)); // Last owned cell
    assert!(!subdomain.is_owned_cell(0, 1, 0)); // Ghost cell
    assert!(!subdomain.is_owned_cell(26, 1, 0)); // Ghost cell

    // Test global index conversion
    let (i_global, j_global, k_global) = subdomain.local_to_global(1, 1, 0);
    assert_eq!(i_global, 25); // i_start_global + (1 - ghost_layers)
    assert_eq!(j_global, 0); // j_start_global + (1 - ghost_layers)
    assert_eq!(k_global, 0); // k_start_global + (0 - ghost_layers)

    // Test ownership checks
    assert!(subdomain.owns_global_cell(25, 0, 0));
    assert!(subdomain.owns_global_cell(49, 19, 0));
    assert!(!subdomain.owns_global_cell(24, 0, 0)); // Previous subdomain
    assert!(!subdomain.owns_global_cell(50, 0, 0)); // Next subdomain
}

#[test]
fn test_domain_decomposition_1d() {
    let global = GlobalExtents::new_2d(100, 50, (0.0, 1.0, 0.0, 1.0));

    // Test with 4 processes
    let subdomain = DomainDecomposition::decompose_1d(
        &global, 1, 4, // rank 1 of 4
    )
    .unwrap();

    assert_eq!(subdomain.rank, 1);
    assert_eq!(subdomain.nx_local, 25); // 100 / 4 = 25 each
    assert_eq!(subdomain.ny_local, 50); // Full y dimension
    assert_eq!(subdomain.i_start_global, 25); // Starts at 25
    assert_eq!(subdomain.ghost_layers, 1);
}

#[test]
fn test_domain_decomposition_2d() {
    let global = GlobalExtents::new_2d(64, 64, (0.0, 1.0, 0.0, 1.0));

    // Test with 4 processes (2x2 grid)
    let subdomain = DomainDecomposition::decompose_2d(
        &global, 2, 4, // rank 2 of 4 (bottom-right in 2x2 grid)
    )
    .unwrap();

    assert_eq!(subdomain.rank, 2);
    assert_eq!(subdomain.nx_local, 32); // 64 / 2 = 32
    assert_eq!(subdomain.ny_local, 32); // 64 / 2 = 32
    assert_eq!(subdomain.i_start_global, 32); // Second column
    assert_eq!(subdomain.j_start_global, 0); // First row
}

#[test]
fn test_utility_functions() {
    // Test that utility functions exist (would need MPI to actually test)
    // These are compile-time checks that the functions exist
    let _ = super::communicator::utils::is_root;
    let _ = super::communicator::utils::barrier;
    let _ = super::communicator::utils::print_root;
    let _ = super::communicator::utils::print_rank;
}

#[test]
fn test_ghost_cell_manager_creation() {
    use std::collections::HashMap;

    // Create mock communicator and neighbors
    let neighbors = HashMap::new();
    // Test that GhostCellManager can be created (compile-time check)
    // In a real test with MPI, we'd need an actual communicator
    let _manager_type: std::marker::PhantomData<GhostCellManager<f64>> = std::marker::PhantomData;
}

#[test]
fn test_neighbor_direction_enum() {
    // Test that all neighbor directions are defined
    let _left = NeighborDirection::Left;
    let _right = NeighborDirection::Right;
    let _bottom = NeighborDirection::Bottom;
    let _top = NeighborDirection::Top;
    let _front = NeighborDirection::Front;
    let _back = NeighborDirection::Back;
}

#[test]
fn test_neighbor_info_creation() {
    let neighbor_info = NeighborInfo {
        direction: NeighborDirection::Left,
        overlap: 1,
    };

    assert_eq!(neighbor_info.overlap, 1);
    match neighbor_info.direction {
        NeighborDirection::Left => {}
        _ => panic!("Wrong direction"),
    }
}

#[test]
fn test_ghost_cell_update_config() {
    let config = GhostCellUpdate {
        max_iterations: 10,
        tolerance: 1e-6,
        blocking: true,
    };

    assert_eq!(config.max_iterations, 10);
    assert_eq!(config.tolerance, 1e-6);
    assert!(config.blocking);
}

#[test]
fn test_ghost_cell_stats_default() {
    let stats = GhostCellStats::default();

    assert_eq!(stats.bytes_sent, 0);
    assert_eq!(stats.bytes_received, 0);
    assert_eq!(stats.messages_sent, 0);
    assert_eq!(stats.messages_received, 0);
    assert_eq!(stats.comm_time, std::time::Duration::default());
}

#[test]
fn test_mpi_error_types() {
    // Test that error types can be created
    let init_error = MpiError::InitializationError("test".to_string());
    let decomp_error = MpiError::DecompositionError("test".to_string());
    let comm_error = MpiError::CommunicationError("test".to_string());
    let not_avail = MpiError::NotAvailable("test".to_string());

    // Check error messages
    match init_error {
        MpiError::InitializationError(msg) => assert_eq!(msg, "test"),
        _ => panic!("Wrong error type"),
    }

    match decomp_error {
        MpiError::DecompositionError(msg) => assert_eq!(msg, "test"),
        _ => panic!("Wrong error type"),
    }

    match comm_error {
        MpiError::CommunicationError(msg) => assert_eq!(msg, "test"),
        _ => panic!("Wrong error type"),
    }

    match not_avail {
        MpiError::NotAvailable(msg) => assert_eq!(msg, "test"),
        _ => panic!("Wrong error type"),
    }
}

#[test]
fn test_domain_decomposition_neighbors_1d() {
    let global = GlobalExtents::new_2d(100, 50, (0.0, 1.0, 0.0, 1.0));
    let subdomain = LocalSubdomain {
        rank: 1,
        nx_local: 25,
        ny_local: 50,
        nz_local: 1,
        i_start_global: 25,
        j_start_global: 0,
        k_start_global: 0,
        ghost_layers: 1,
    };

    // Test neighbor computation for middle subdomain in 1D decomposition
    let neighbors = DomainDecomposition::compute_neighbors(
        &subdomain,
        &global,
        4,
        DecompositionStrategy::Simple1D,
    )
    .unwrap();

    // Should have both left and right neighbors
    assert_eq!(neighbors.len(), 2);
    assert!(neighbors.contains_key(&0)); // left neighbor (rank 0)
    assert!(neighbors.contains_key(&2)); // right neighbor (rank 2)

    if let Some(left_info) = neighbors.get(&0) {
        assert_eq!(left_info.direction, NeighborDirection::Left);
        assert_eq!(left_info.overlap, 1);
    }

    if let Some(right_info) = neighbors.get(&2) {
        assert_eq!(right_info.direction, NeighborDirection::Right);
        assert_eq!(right_info.overlap, 1);
    }
}

#[test]
fn test_domain_decomposition_neighbors_boundary() {
    let global = GlobalExtents::new_2d(100, 50, (0.0, 1.0, 0.0, 1.0));
    let subdomain = LocalSubdomain {
        rank: 0, // Leftmost subdomain
        nx_local: 25,
        ny_local: 50,
        nz_local: 1,
        i_start_global: 0,
        j_start_global: 0,
        k_start_global: 0,
        ghost_layers: 1,
    };

    // Test neighbor computation for left boundary subdomain
    let neighbors = DomainDecomposition::compute_neighbors(
        &subdomain,
        &global,
        4,
        DecompositionStrategy::Simple1D,
    )
    .unwrap();

    // Should only have right neighbor (no left neighbor for boundary)
    assert_eq!(neighbors.len(), 1);
    assert!(neighbors.contains_key(&1)); // right neighbor (rank 1)

    if let Some(right_info) = neighbors.get(&1) {
        assert_eq!(right_info.direction, NeighborDirection::Right);
        assert_eq!(right_info.overlap, 1);
    }
}

#[test]
fn test_find_rank_for_global_index() {
    let global = GlobalExtents::new_2d(100, 50, (0.0, 1.0, 0.0, 1.0));

    // Test rank finding for 1D decomposition with 4 processes
    // Each process gets 25 cells: [0-24], [25-49], [50-74], [75-99]

    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(10, &global, 4).unwrap(),
        0
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(30, &global, 4).unwrap(),
        1
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(60, &global, 4).unwrap(),
        2
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(85, &global, 4).unwrap(),
        3
    );

    // Test boundary cases
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(24, &global, 4).unwrap(),
        0
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(25, &global, 4).unwrap(),
        1
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(49, &global, 4).unwrap(),
        1
    );
    assert_eq!(
        DomainDecomposition::find_rank_for_global_i(50, &global, 4).unwrap(),
        2
    );
}

#[test]
fn test_distributed_linear_operator_trait() {
    // Compile-time test that trait can be implemented

    // Test that the trait exists and has required methods
    let _local_dim: usize = 100;
    let _global_dim: usize = 1000;

    // Test trait bounds
    fn _test_trait<T: RealField, Op: DistributedLinearOperator<T>>(_op: &Op) {
        let _local = _op.local_dimension();
        let _global = _op.global_dimension();
    }
}

#[test]
fn test_parallel_preconditioner_types() {
    // Test that preconditioner types can be instantiated (compile-time check)

    let _jacobi_type: std::marker::PhantomData<
        BlockJacobiPreconditioner<f64, DistributedLaplacian2D<f64>>,
    > = std::marker::PhantomData;

    let _schwarz_type: std::marker::PhantomData<
        AdditiveSchwarzPreconditioner<f64, DistributedLaplacian2D<f64>>,
    > = std::marker::PhantomData;
}

#[test]
fn test_distributed_gmres_solver_type() {
    // Test that DistributedGMRES can be instantiated (compile-time check)

    let _gmres_type: std::marker::PhantomData<
        DistributedGMRES<
            f64,
            DistributedLaplacian2D<f64>,
            BlockJacobiPreconditioner<f64, DistributedLaplacian2D<f64>>,
        >,
    > = std::marker::PhantomData;
}

#[test]
fn test_parallel_io_types() {
    // Test that parallel I/O types exist and can be instantiated
    use super::distributed_solvers::parallel_io::ParallelVtkWriter;

    let _vtk_writer: std::marker::PhantomData<ParallelVtkWriter<f64>> = std::marker::PhantomData;
}

#[test]
fn test_load_balance_metrics_creation() {
    let metrics = LoadBalanceMetrics {
        max_load: 1000,
        min_load: 800,
        avg_load: 900.0,
        imbalance_ratio: 1.11,
        efficiency: 0.9,
        needs_rebalancing: false,
    };

    assert_eq!(metrics.max_load, 1000);
    assert_eq!(metrics.min_load, 800);
    assert_eq!(metrics.avg_load, 900.0);
    assert_eq!(metrics.imbalance_ratio, 1.11);
    assert_eq!(metrics.efficiency, 0.9);
    assert!(!metrics.needs_rebalancing);
}

#[test]
fn test_refinement_criteria_creation() {
    let criteria = RefinementCriteria {
        error_threshold: 1e-3,
        coarsening_threshold: 1e-5,
        max_refinement_ratio: 4,
    };

    assert_eq!(criteria.error_threshold, 1e-3);
    assert_eq!(criteria.coarsening_threshold, 1e-5);
    assert_eq!(criteria.max_refinement_ratio, 4);
}

#[test]
fn test_adaptive_mesh_refinement_creation() {
    let criteria = RefinementCriteria {
        error_threshold: 1e-3,
        coarsening_threshold: 1e-5,
        max_refinement_ratio: 4,
    };

    let amr = AdaptiveMeshRefinement::new(5, criteria, None);

    assert_eq!(amr.refinement_level(), 0);
    assert_eq!(amr.max_refinement_level, 5);
}

#[test]
fn test_adaptive_mesh_refinement_needs_refinement() {
    let criteria = RefinementCriteria {
        error_threshold: 1e-3,
        coarsening_threshold: 1e-5,
        max_refinement_ratio: 4,
    };

    let amr = AdaptiveMeshRefinement::new(5, criteria, None);

    let error_estimates = vec![1e-4, 2e-3, 5e-4, 1e-2]; // Errors: below, above, below, above threshold
    let needs_refine = amr.needs_refinement(&error_estimates);

    assert_eq!(needs_refine, vec![false, true, false, true]);
}

#[test]
fn test_adaptive_mesh_refinement_needs_coarsening() {
    let criteria = RefinementCriteria {
        error_threshold: 1e-3,
        coarsening_threshold: 1e-5,
        max_refinement_ratio: 4,
    };

    let mut amr = AdaptiveMeshRefinement::new(5, criteria, None);
    amr.increment_level(); // Set level to 1 to enable coarsening

    let error_estimates = vec![1e-6, 1e-4, 1e-7, 1e-3]; // Errors: below, above, below, above threshold
    let needs_coarsen = amr.needs_coarsening(&error_estimates);

    assert_eq!(needs_coarsen, vec![true, false, true, false]);
}

#[test]
fn test_adaptive_mesh_refinement_level_management() {
    let criteria = RefinementCriteria {
        error_threshold: 1e-3,
        coarsening_threshold: 1e-5,
        max_refinement_ratio: 4,
    };

    let mut amr = AdaptiveMeshRefinement::new(3, criteria, None);

    assert_eq!(amr.refinement_level(), 0);

    amr.increment_level();
    assert_eq!(amr.refinement_level(), 1);

    amr.increment_level();
    assert_eq!(amr.refinement_level(), 2);

    amr.increment_level();
    assert_eq!(amr.refinement_level(), 3);

    // Should not exceed max level
    amr.increment_level();
    assert_eq!(amr.refinement_level(), 3);
}

#[test]
fn test_load_balance_metrics_calculation() {
    // Test load balance metrics with perfect balance
    let metrics = LoadBalanceMetrics {
        max_load: 1000,
        min_load: 1000,
        avg_load: 1000.0,
        imbalance_ratio: 1.0,
        efficiency: 1.0,
        needs_rebalancing: false,
    };

    assert_eq!(metrics.imbalance_ratio, 1.0);
    assert_eq!(metrics.efficiency, 1.0);
    assert!(!metrics.needs_rebalancing);

    // Test with imbalance
    let imbalanced_metrics = LoadBalanceMetrics {
        max_load: 1500,
        min_load: 500,
        avg_load: 1000.0,
        imbalance_ratio: 1.5,
        efficiency: 2.0 / 3.0,
        needs_rebalancing: true,
    };

    assert_eq!(imbalanced_metrics.imbalance_ratio, 1.5);
    assert!((imbalanced_metrics.efficiency - 2.0 / 3.0).abs() < 1e-10);
    assert!(imbalanced_metrics.needs_rebalancing);
}

#[test]
fn test_distributed_vector_creation() {
    // This would require MPI initialization in a real test
    // Placeholder for compile-time validation
    let _marker: std::marker::PhantomData<DistributedVector<f64>> = std::marker::PhantomData;
}
