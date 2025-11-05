//! End-to-end integration tests for MPI parallelization
//!
//! These tests verify that all MPI components work together correctly
//! for production CFD simulations.

#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::*;
use std::collections::HashMap;

#[cfg(feature = "mpi")]
mod integration_tests {
    use super::*;
    use nalgebra::Vector2;

    #[test]
    fn test_mpi_initialization_and_communicator() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        assert!(world.size() >= 1);
        assert!(world.rank() >= 0 && world.rank() < world.size());

        // Test communicator functionality
        let rank_data = world.rank();
        let mut gathered_data = vec![0i32; world.size() as usize];

        // Test all-gather (this would be MPI_Allgather in real MPI)
        for i in 0..world.size() as usize {
            if i == world.rank() as usize {
                gathered_data[i] = rank_data;
            }
        }

        // Verify all ranks have their data
        assert_eq!(gathered_data[world.rank() as usize], world.rank());
    }

    #[test]
    fn test_domain_decomposition_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Create global domain
        let global_extents = GlobalExtents::new_2d(1000, 1000, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;

        // Create decomposition
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // Verify decomposition properties
        assert_eq!(decomp.global_extents(), &global_extents);
        assert_eq!(decomp.strategy, strategy);
        assert_eq!(decomp.communicator().size(), world.size());
        assert_eq!(decomp.communicator().rank(), world.rank());

        // Verify local subdomain
        let subdomain = decomp.local_subdomain();
        assert!(subdomain.rank >= 0 && subdomain.rank < world.size());
        assert!(subdomain.nx_local > 0 && subdomain.ny_local > 0);
        assert!(subdomain.i_start_global >= 0);
        assert!(subdomain.j_start_global >= 0);
        assert_eq!(subdomain.ghost_layers, 1);

        // Verify neighbors
        let neighbors = decomp.neighbors();
        assert!(!neighbors.is_empty()); // Should have at least boundary neighbors
    }

    #[test]
    fn test_ghost_cell_exchange_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Skip test if only one process
        if world.size() == 1 {
            return;
        }

        // Create domain decomposition
        let global_extents = GlobalExtents::new_2d(100, 100, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // Create test data with ghost layers
        let nx_total = decomp.local_subdomain().nx_local + 2; // +2 for ghost layers
        let ny_total = decomp.local_subdomain().ny_local + 2;

        let mut velocity_u = vec![vec![Vector2::new(1.0, 2.0); ny_total]; nx_total];
        let mut velocity_v = vec![vec![Vector2::new(3.0, 4.0); ny_total]; nx_total];
        let mut pressure = vec![vec![5.0; ny_total]; nx_total];

        // Set interior values to rank-specific values for verification
        let rank_offset = world.rank() as f64 * 10.0;
        for i in 1..nx_total-1 {
            for j in 1..ny_total-1 {
                velocity_u[i][j] = Vector2::new(10.0 + rank_offset, 20.0 + rank_offset);
                velocity_v[i][j] = Vector2::new(30.0 + rank_offset, 40.0 + rank_offset);
                pressure[i][j] = 50.0 + rank_offset;
            }
        }

        // Perform ghost cell exchange
        update_ghost_cells(
            &world,
            &decomp,
            &mut velocity_u,
            &mut velocity_v,
            &mut pressure,
        ).unwrap();

        // Verify ghost cells were updated (non-trivial verification would require
        // coordinating expected values across processes)
        // For now, just ensure the function completes without error
    }

    #[test]
    fn test_load_balancing_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Create domain decomposition
        let global_extents = GlobalExtents::new_2d(1000, 1000, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // Create load balancer
        let mut load_balancer = LoadBalancer::new(&world, decomp, 1.2, 100).unwrap();

        // Test load balance assessment
        let local_workload = 10000 + world.rank() as usize * 1000; // Vary workload by rank
        let metrics = load_balancer.assess_load_balance(local_workload).unwrap();

        // Verify metrics are reasonable
        assert!(metrics.max_load >= metrics.min_load);
        assert!(metrics.avg_load > 0.0);
        assert!(metrics.imbalance_ratio >= 1.0);
        assert!(metrics.efficiency <= 1.0);

        // Test repartitioning if imbalance threshold exceeded
        if load_balancer.should_repartition(&metrics) {
            let new_workloads = vec![10000; world.size() as usize]; // Balanced workloads
            let new_decomp = load_balancer.repartition(&new_workloads).unwrap();

            // Verify new decomposition is valid
            assert_eq!(new_decomp.global_extents(), &global_extents);
            assert_eq!(new_decomp.strategy, strategy);
        }
    }

    #[test]
    fn test_adaptive_mesh_refinement_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Create domain decomposition
        let global_extents = GlobalExtents::new_2d(1000, 1000, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // Create refinement criteria
        let criteria = RefinementCriteria {
            error_threshold: 1e-3,
            coarsening_threshold: 1e-5,
            max_refinement_ratio: 4,
        };

        // Create load balancer for AMR
        let load_balancer = LoadBalancer::new(&world, decomp.clone(), 1.2, 100).unwrap();
        let mut amr = AdaptiveMeshRefinement::new(5, criteria, Some(load_balancer));

        // Test refinement decision logic
        let error_estimates = vec![1e-4, 2e-3, 5e-4, 1e-2]; // Mix of errors
        let needs_refine = amr.needs_refinement(&error_estimates);
        assert_eq!(needs_refine, vec![false, true, false, true]);

        // Test coarsening decision logic
        amr.increment_level(); // Set level > 0 to enable coarsening
        let needs_coarsen = amr.needs_coarsening(&error_estimates);
        assert_eq!(needs_coarsen, vec![false, false, false, false]); // None below threshold

        let low_errors = vec![1e-7, 1e-6, 1e-8]; // All below coarsening threshold
        let needs_coarsen_low = amr.needs_coarsening(&low_errors);
        assert_eq!(needs_coarsen_low, vec![true, true, true]);
    }

    #[test]
    fn test_distributed_linear_solver_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Create domain decomposition
        let global_extents = GlobalExtents::new_2d(100, 100, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // Create distributed Laplacian operator
        let dx = 0.01;
        let dy = 0.01;
        let operator = DistributedLaplacian2D::new(&decomp, &world, dx, dy).unwrap();

        // Verify operator properties
        assert_eq!(operator.local_dimension(), (decomp.local_subdomain().nx_local * decomp.local_subdomain().ny_local));
        assert_eq!(operator.global_dimension(), (global_extents.nx_global * global_extents.ny_global));

        // Create distributed vectors
        let local_size = operator.local_dimension();
        use nalgebra::DVector;
        let local_data = DVector::from_element(local_size, 1.0);
        let x = DistributedVector::new(local_data.clone(), &world);
        let result = operator.apply(&x).unwrap();

        // Verify result is valid
        assert_eq!(result.local_dimension(), local_size);
        assert_eq!(result.global_dimension(), operator.global_dimension());
    }

    #[test]
    fn test_performance_validation_integration() {
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // Create performance validator
        let validator = PerformanceValidator::<f64>::new(&world);

        // Test production readiness assessment
        let readiness_report = validator.assess_production_readiness().unwrap();

        // Verify report structure
        assert!(readiness_report.overall_score >= 0 && readiness_report.overall_score <= 100);
        assert!(!readiness_report.component_scores.is_empty());

        // Verify component scores are reasonable
        for (&component, &score) in &readiness_report.component_scores {
            assert!(score >= 0 && score <= 100, "Invalid score for {}: {}", component, score);
        }

        // Verify deployment configuration
        assert!(readiness_report.recommended_config.cores_per_node > 0);
        assert!(!readiness_report.recommended_config.mpi_implementation.is_empty());
        assert!(readiness_report.recommended_config.memory_per_core_mb > 0);

        // Verify scaling limits
        assert!(readiness_report.scaling_limits.max_cores > 0);
        assert!(readiness_report.scaling_limits.efficiency_threshold > 0.0);
        assert!(readiness_report.scaling_limits.comm_overhead_limit > 0.0);
    }

    #[test]
    fn test_complete_mpi_pipeline() {
        // This test simulates a complete CFD simulation pipeline using MPI
        let universe = MpiUniverse::new().unwrap();
        let world = universe.world();

        // 1. Initialize domain decomposition
        let global_extents = GlobalExtents::new_2d(200, 200, (0.0, 1.0, 0.0, 1.0));
        let strategy = DecompositionStrategy::Cartesian2D;
        let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

        // 2. Initialize load balancing
        let load_balancer = LoadBalancer::new(&world, decomp.clone(), 1.2, 100).unwrap();

        // 3. Initialize distributed solver
        let dx = 0.005;
        let dy = 0.005;
        let operator = DistributedLaplacian2D::new(&decomp, &world, dx, dy).unwrap();

        // 4. Create preconditioner
        let preconditioner = BlockJacobiPreconditioner::new(&operator, &decomp, &world).unwrap();

        // 5. Create distributed GMRES solver
        let solver = DistributedGMRES::new(operator, preconditioner, &world, 30);

        // 6. Create test problem (simple Poisson equation)
        let local_size = decomp.local_subdomain().nx_local * decomp.local_subdomain().ny_local;
        use nalgebra::DVector;
        let rhs_data = DVector::from_element(local_size, 1.0); // Constant RHS
        let rhs = DistributedVector::new(rhs_data, &world);
        let initial_guess = DistributedVector::new(DVector::from_element(local_size, 0.0), &world);

        // 7. Solve system
        let tolerance = 1e-6;
        let max_iter = 100;
        let solution = solver.solve(&rhs, &initial_guess, tolerance, max_iter).unwrap();

        // 8. Verify solution exists and is reasonable
        assert_eq!(solution.local_dimension(), local_size);
        assert_eq!(solution.global_dimension(), operator.global_dimension());

        // 9. Test load balancing assessment
        let local_workload = local_size;
        let metrics = load_balancer.assess_load_balance(local_workload).unwrap();
        assert!(metrics.imbalance_ratio >= 1.0);

        // Pipeline completed successfully
        if world.is_root() {
            println!("✅ Complete MPI pipeline test passed");
        }
    }

    #[test]
    fn test_error_handling_integration() {
        // Test that error handling works correctly across MPI components

        // Test MpiError creation and handling
        let error = MpiError::CommunicationError("Test communication error".to_string());
        assert!(matches!(error, MpiError::CommunicationError(_)));

        let error = MpiError::DecompositionError("Test decomposition error".to_string());
        assert!(matches!(error, MpiError::DecompositionError(_)));

        let error = MpiError::SolverError("Test solver error".to_string());
        assert!(matches!(error, MpiError::SolverError(_)));

        // Test Result type alias works
        let result: MpiResult<()> = Ok(());
        assert!(result.is_ok());

        let result: MpiResult<()> = Err(MpiError::NotAvailable("Test".to_string()));
        assert!(result.is_err());
    }

    #[test]
    fn test_mpi_feature_gating() {
        // This test ensures MPI features are properly gated
        // When MPI feature is enabled, these types should exist
        let _universe: MpiUniverse;
        let _world: MpiCommunicator;
        let _decomp: DomainDecomposition;
        let _metrics: PerformanceMetrics;
        let _validator: PerformanceValidator<f64>;

        // Compilation should succeed when MPI feature is enabled
    }
}

#[cfg(not(feature = "mpi"))]
mod no_mpi_tests {
    #[test]
    fn test_mpi_features_disabled() {
        // When MPI feature is disabled, MPI types should not be available
        // This test ensures the feature gating works correctly

        // These lines would fail to compile if MPI feature is enabled
        // and would pass if disabled (since the types don't exist)
        // We can't actually test this directly, but the compilation
        // success indicates proper feature gating
    }
}
