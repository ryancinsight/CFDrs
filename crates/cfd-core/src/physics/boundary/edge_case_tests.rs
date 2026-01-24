//! Comprehensive edge case tests for boundary conditions
//!
//! This module provides rigorous testing of boundary condition applicators
//! with property-based testing (proptest) and edge case validation per
//! CFD standards (Patankar 1980, Versteeg & Malalasekera 2007).

#[cfg(test)]
mod boundary_edge_cases {
    use crate::physics::boundary::{
        BoundaryCondition, BoundaryConditionApplicator, BoundaryConditionSpec, DirichletApplicator,
        NeumannApplicator, RobinApplicator,
    };
    use proptest::prelude::*;

    /// Test Dirichlet boundary condition with zero values
    /// Edge case: Zero boundary values should be handled correctly
    #[test]
    fn test_dirichlet_zero_value() {
        let mut field = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 0.0,
                component_values: None,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        assert_eq!(field[0], 0.0, "West boundary should be zero");
        assert_eq!(field[1], 2.0, "Interior should be unchanged");
    }

    /// Test Dirichlet with very large values
    /// Edge case: Large values (O(10^6)) should not cause overflow
    #[test]
    fn test_dirichlet_large_values() {
        let mut field = vec![1.0f64; 10];
        let applicator = DirichletApplicator::<f64>::new();

        let large_value = 1e6;
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: large_value,
                component_values: None,
            },
            "east".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        assert_eq!(
            field[9], large_value,
            "Large value should be applied correctly"
        );
        assert!(field[9].is_finite(), "Value should remain finite");
    }

    /// Test Dirichlet with very small values
    /// Edge case: Small values (O(10^-10)) should maintain precision
    #[test]
    fn test_dirichlet_small_values() {
        let mut field = vec![1.0f64; 10];
        let applicator = DirichletApplicator::<f64>::new();

        let small_value = 1e-10;
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: small_value,
                component_values: None,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        assert!(
            (field[0] - small_value).abs() < 1e-15,
            "Small value precision should be maintained"
        );
    }

    /// Test Dirichlet with negative values
    /// Edge case: Negative boundary values (e.g., reverse flow)
    #[test]
    fn test_dirichlet_negative_values() {
        let mut field = vec![1.0f64; 5];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: -10.0,
                component_values: None,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        assert_eq!(field[0], -10.0, "Negative values should be applied");
    }

    /// Test Dirichlet on empty field
    /// Edge case: Empty field should not panic
    #[test]
    fn test_dirichlet_empty_field() {
        let mut field: Vec<f64> = vec![];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 5.0,
                component_values: None,
            },
            "west".to_string(),
        );

        let result = applicator.apply(&mut field, &spec, 0.0);
        assert!(result.is_ok(), "Empty field should not cause error");
    }

    /// Test Dirichlet on single-element field
    /// Edge case: Single element field (both boundaries are same point)
    #[test]
    fn test_dirichlet_single_element() {
        let mut field = vec![1.0f64];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 10.0,
                component_values: None,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        assert_eq!(field[0], 10.0, "Single element should be updated");
    }

    /// Test Neumann boundary condition with zero gradient
    /// Edge case: Zero gradient (insulated boundary)
    #[test]
    fn test_neumann_zero_gradient() {
        let mut field = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let applicator = NeumannApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Neumann { gradient: 0.0 },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        // Zero gradient means boundary = interior neighbor
        assert_eq!(
            field[0], field[1],
            "Zero gradient: boundary equals neighbor"
        );
    }

    /// Test Neumann with large positive gradient
    /// Edge case: Large gradient values
    #[test]
    fn test_neumann_large_gradient() {
        let mut field = vec![1.0f64; 10];
        let applicator = NeumannApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Neumann { gradient: 1000.0 },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        // Neumann: u_boundary = u_interior - dx * gradient
        // With dx ~ 1.0, boundary should be significantly different
        assert_ne!(field[0], field[1], "Large gradient should affect boundary");
    }

    /// Test Neumann with negative gradient
    /// Edge case: Negative gradient (heat flux out)
    #[test]
    fn test_neumann_negative_gradient() {
        let mut field = vec![10.0f64; 5];
        let applicator = NeumannApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Neumann { gradient: -5.0 },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        // Negative gradient should increase boundary value
        assert!(field[0] >= field[1], "Negative gradient effect");
    }

    /// Test Robin boundary condition with various mixing ratios
    /// Edge case: Robin BC interpolates between Dirichlet and Neumann
    #[test]
    fn test_robin_mixing_ratios() {
        // Test α=1 (pure Dirichlet-like)
        let mut field = vec![1.0f64, 2.0, 3.0];
        let applicator = RobinApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Robin {
                alpha: 1.0,
                beta: 0.0,
                gamma: 5.0,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();
        // α*u + β*du/dx = γ with β=0 => u = γ/α = 5.0
        assert!(
            (field[0] - 5.0).abs() < 1e-10,
            "Robin α=1 should behave like Dirichlet"
        );
    }

    /// Test Robin with zero alpha (pure Neumann-like)
    /// Edge case: α→0 degenerates to Neumann
    #[test]
    fn test_robin_zero_alpha() {
        let mut field = vec![1.0f64, 2.0, 3.0];
        let applicator = RobinApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Robin {
                alpha: 1e-10, // Nearly zero
                beta: 1.0,
                gamma: 0.0,
            },
            "west".to_string(),
        );

        let result = applicator.apply(&mut field, &spec, 0.0);
        assert!(result.is_ok(), "Small alpha should not cause issues");
    }

    // Property-based test: Dirichlet preserves value exactly
    // Property: Applied value should equal specified value (within FP precision)
    proptest! {
        #[test]
        fn prop_dirichlet_preserves_value(value in -1000.0..1000.0f64) {
            let mut field = vec![0.0f64; 10];
            let applicator = DirichletApplicator::<f64>::new();

            let spec = BoundaryConditionSpec::new(
                BoundaryCondition::Dirichlet {
                    value,
                    component_values: None,
                },
                "west".to_string(),
            );

            applicator.apply(&mut field, &spec, 0.0).unwrap();

            // Property: Boundary value equals specified value
            assert!((field[0] - value).abs() < 1e-10, "Dirichlet must preserve value exactly");
        }
    }

    // Property-based test: Interior points unchanged by boundary application
    // Property: Boundary conditions should not modify interior points
    proptest! {
        #[test]
        fn prop_interior_unchanged(
            value in -100.0..100.0f64,
            interior_val in -100.0..100.0f64
        ) {
            let mut field = vec![interior_val; 10];
            let applicator = DirichletApplicator::<f64>::new();

            let spec = BoundaryConditionSpec::new(
                BoundaryCondition::Dirichlet {
                    value,
                    component_values: None,
                },
                "west".to_string(),
            );

            applicator.apply(&mut field, &spec, 0.0).unwrap();

            // Property: Interior points (indices 1..9) should be unchanged
            for &val in field.iter().take(9).skip(1) {
                assert_eq!(val, interior_val, "Interior should remain unchanged");
            }
        }
    }

    // Property-based test: Neumann zero gradient creates constant field at boundary
    // Property: Zero gradient => boundary value equals neighbor
    proptest! {
        #[test]
        fn prop_neumann_zero_gradient_constant(field_val in -100.0..100.0f64) {
            let mut field = vec![field_val; 10];
            let applicator = NeumannApplicator::<f64>::new();

            let spec = BoundaryConditionSpec::new(
                BoundaryCondition::Neumann { gradient: 0.0 },
                "west".to_string(),
            );

            applicator.apply(&mut field, &spec, 0.0).unwrap();

            // Property: Zero gradient means field[0] == field[1]
            assert_eq!(field[0], field[1], "Zero gradient should create continuity");
        }
    }

    // Property-based test: Robin BC interpolation property
    // Property: α=1, β=0 behaves like Dirichlet
    proptest! {
        #[test]
        fn prop_robin_dirichlet_equivalence(gamma in -100.0..100.0f64) {
            let mut field = vec![1.0f64; 5];
            let applicator = RobinApplicator::<f64>::new();

            let spec = BoundaryConditionSpec::new(
                BoundaryCondition::Robin {
                    alpha: 1.0,
                    beta: 0.0,
                    gamma,
                },
                "west".to_string(),
            );

            applicator.apply(&mut field, &spec, 0.0).unwrap();

            // Property: α*u = γ with α=1 => u = γ
            assert!((field[0] - gamma).abs() < 1e-10, "Robin should reduce to Dirichlet");
        }
    }

    /// Edge case: Multiple boundary applications should be idempotent
    /// Property: Applying same BC twice yields same result
    #[test]
    fn test_boundary_idempotence() {
        let mut field1 = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let mut field2 = field1.clone();
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 10.0,
                component_values: None,
            },
            "west".to_string(),
        );

        // Apply once
        applicator.apply(&mut field1, &spec, 0.0).unwrap();

        // Apply twice
        applicator.apply(&mut field2, &spec, 0.0).unwrap();
        applicator.apply(&mut field2, &spec, 0.0).unwrap();

        // Should yield same result
        assert_eq!(field1, field2, "Boundary application should be idempotent");
    }

    /// Edge case: Time-dependent boundary condition evaluation
    /// Boundary value should be constant for constant BC
    #[test]
    fn test_constant_time_evaluation() {
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 5.0,
                component_values: None,
            },
            "west".to_string(),
        );

        // Constant BC should be same at all times
        let bc_t0 = spec.evaluate_at_time(0.0);
        let bc_t1 = spec.evaluate_at_time(1.0);

        assert_eq!(
            bc_t0.as_ref(),
            bc_t1.as_ref(),
            "Constant BC should not vary with time"
        );
    }

    /// Edge case: Very large field sizes
    /// Test boundary application on large fields (performance check)
    #[test]
    fn test_large_field_boundary() {
        let size = 10_000;
        let mut field = vec![1.0f64; size];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 42.0,
                component_values: None,
            },
            "west".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();

        assert_eq!(field[0], 42.0, "Boundary should be applied on large field");
        assert_eq!(field[size - 1], 1.0, "Interior should be unchanged");
    }

    /// Edge case: All-boundary region application
    /// Test applying condition to all elements
    #[test]
    fn test_all_boundary_application() {
        let mut field = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 10.0,
                component_values: None,
            },
            "all".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();

        // All elements should be set to 10.0
        for val in &field {
            assert_eq!(*val, 10.0, "All elements should be updated");
        }
    }

    /// Edge case: East boundary application
    /// Test applying condition to east (right) boundary
    #[test]
    fn test_east_boundary_application() {
        let mut field = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let applicator = DirichletApplicator::<f64>::new();

        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 20.0,
                component_values: None,
            },
            "east".to_string(),
        );

        applicator.apply(&mut field, &spec, 0.0).unwrap();

        assert_eq!(field[0], 1.0, "West boundary should be unchanged");
        assert_eq!(field[4], 20.0, "East boundary should be updated");
    }

    /// Edge case: Condition type identification
    /// Test that condition types are correctly identified
    #[test]
    fn test_condition_type_identification() {
        let dirichlet_spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet {
                value: 5.0,
                component_values: None,
            },
            "west".to_string(),
        );
        assert_eq!(dirichlet_spec.condition_type(), "Dirichlet");

        let neumann_spec = BoundaryConditionSpec::new(
            BoundaryCondition::Neumann { gradient: 1.0 },
            "east".to_string(),
        );
        assert_eq!(neumann_spec.condition_type(), "Neumann");

        let robin_spec = BoundaryConditionSpec::new(
            BoundaryCondition::Robin {
                alpha: 1.0,
                beta: 1.0,
                gamma: 0.0,
            },
            "north".to_string(),
        );
        assert_eq!(robin_spec.condition_type(), "Robin");
    }

    /// Edge case: Applicator supports check
    /// Test that applicators correctly identify supported conditions
    #[test]
    fn test_applicator_supports() {
        let dirichlet_condition = BoundaryCondition::Dirichlet {
            value: 5.0,
            component_values: None,
        };
        let neumann_condition = BoundaryCondition::Neumann { gradient: 1.0 };

        let dirichlet_applicator = DirichletApplicator::<f64>::new();
        let neumann_applicator = NeumannApplicator::<f64>::new();

        assert!(
            dirichlet_applicator.supports(&dirichlet_condition),
            "Dirichlet applicator should support Dirichlet condition"
        );
        assert!(
            !dirichlet_applicator.supports(&neumann_condition),
            "Dirichlet applicator should not support Neumann condition"
        );
        assert!(
            neumann_applicator.supports(&neumann_condition),
            "Neumann applicator should support Neumann condition"
        );
        assert!(
            !neumann_applicator.supports(&dirichlet_condition),
            "Neumann applicator should not support Dirichlet condition"
        );
    }

    /// Edge case: Applicator name check
    /// Test that applicators report correct names
    #[test]
    fn test_applicator_names() {
        let dirichlet_applicator = DirichletApplicator::<f64>::new();
        let neumann_applicator = NeumannApplicator::<f64>::new();
        let robin_applicator = RobinApplicator::<f64>::new();

        assert_eq!(dirichlet_applicator.name(), "Dirichlet");
        assert_eq!(neumann_applicator.name(), "Neumann");
        assert_eq!(robin_applicator.name(), "Robin");
    }
}
