//! Hierarchical integration harness for validation contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (37 binaries linking the stack per build become 3: this harness plus
//! the two allocator singletons below). Nextest still isolates per test,
//! and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.
//!
//! `allocator_compat` and `tracking_allocator` stay standalone binaries by
//! contract: each installs its own process-global allocator (`System` vs
//! the tracking instrument), and one binary admits exactly one
//! `#[global_allocator]`. Their headers document the solitude requirement.

#[path = "main/advanced_physics_validation.rs"]
mod advanced_physics_validation;
#[path = "main/analytical_formulas.rs"]
mod analytical_formulas;
#[path = "main/analytical_poiseuille_3d.rs"]
mod analytical_poiseuille_3d;
#[path = "main/automated_validation_suite.rs"]
mod automated_validation_suite;
#[path = "main/benchmark_validation.rs"]
mod benchmark_validation;
#[path = "main/cell_separation_validation.rs"]
mod cell_separation_validation;
#[path = "main/complete_validation_suite.rs"]
mod complete_validation_suite;
#[path = "main/complex_boundary_mms_validation.rs"]
mod complex_boundary_mms_validation;
#[path = "main/comprehensive_validation_pipeline.rs"]
mod comprehensive_validation_pipeline;
#[path = "main/cross_fidelity_1d_2d_3d.rs"]
mod cross_fidelity_1d_2d_3d;
#[path = "main/cross_fidelity_1d_vs_2d.rs"]
mod cross_fidelity_1d_vs_2d;
#[path = "main/cross_fidelity_branching_public_api.rs"]
mod cross_fidelity_branching_public_api;
#[path = "main/cross_fidelity_circular_duct.rs"]
mod cross_fidelity_circular_duct;
#[path = "main/cross_fidelity_non_newtonian.rs"]
mod cross_fidelity_non_newtonian;
#[path = "main/cross_fidelity_serpentine.rs"]
mod cross_fidelity_serpentine;
#[path = "main/cross_fidelity_shape_factors.rs"]
mod cross_fidelity_shape_factors;
#[path = "main/cross_fidelity_trifurcation.rs"]
mod cross_fidelity_trifurcation;
#[path = "main/cross_fidelity_venturi_calibration.rs"]
mod cross_fidelity_venturi_calibration;
#[path = "main/enhanced_benchmarking.rs"]
mod enhanced_benchmarking;
#[path = "main/integration_tests.rs"]
mod integration_tests;
#[path = "main/literature_benchmarks.rs"]
mod literature_benchmarks;
#[path = "main/mms_comprehensive_validation.rs"]
mod mms_comprehensive_validation;
#[path = "main/mms_edge_cases.rs"]
mod mms_edge_cases;
#[path = "main/multi_layer_junction_validation.rs"]
mod multi_layer_junction_validation;
#[path = "main/multi_physics_validation.rs"]
mod multi_physics_validation;
#[path = "main/physics_model_validation.rs"]
mod physics_model_validation;
#[path = "main/physics_validation.rs"]
mod physics_validation;
#[path = "main/proptest_convergence.rs"]
mod proptest_convergence;
#[path = "main/richardson_third_order.rs"]
mod richardson_third_order;
#[path = "main/schematics2mesh.rs"]
mod schematics2mesh;
#[path = "main/therapy_bifurcation_val.rs"]
mod therapy_bifurcation_val;
#[path = "main/turbulence_model_validation.rs"]
mod turbulence_model_validation;
#[path = "main/turbulent_mms_validation.rs"]
mod turbulent_mms_validation;
#[path = "main/twelve_steps.rs"]
mod twelve_steps;
#[path = "main/vortex_shear_validation.rs"]
mod vortex_shear_validation;
