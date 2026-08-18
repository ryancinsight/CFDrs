#![allow(clippy::float_cmp)]
#![allow(clippy::print_stdout)]
#![allow(clippy::print_stderr)]
use super::{AdvectionConfig, GpuAdvectionKernel};
use crate::compute::gpu::GpuContext;
use crate::error::Error;
use std::sync::Arc;

fn kernel() -> Option<GpuAdvectionKernel> {
    let context = match GpuContext::create() {
        Ok(context) => context,
        Err(error) => {
            eprintln!("Skipping GPU advection test: {error:?}");
            return None;
        }
    };
    match GpuAdvectionKernel::new(Arc::new(context)) {
        Ok(kernel) => Some(kernel),
        Err(error) => {
            eprintln!("Skipping GPU advection test: {error:?}");
            None
        }
    }
}

#[test]
fn zero_velocity_is_exact_identity_across_partial_workgroups_and_planes() {
    let Some(kernel) = kernel() else { return };
    let dimensions = [9, 5, 2];
    let config = AdvectionConfig::new(dimensions, [1.0; 3], 1.0).expect("expected value");
    let scalar: Vec<f32> = (0..config.element_count())
        .map(|index| index as f32)
        .collect();
    let velocity = vec![0.0; config.element_count()];
    let mut output = vec![-1.0; config.element_count()];

    kernel
        .execute(&scalar, &velocity, &velocity, config, &mut output)
        .expect("expected value");

    assert_eq!(output, scalar);
}

#[test]
fn directional_upwind_selection_is_exact_and_boundaries_are_copied() {
    let Some(kernel) = kernel() else { return };
    let dimensions = [5, 4, 1];
    let config = AdvectionConfig::new(dimensions, [1.0; 3], 0.25).expect("expected value");
    let scalar: Vec<f32> = (0..dimensions[1])
        .flat_map(|_| (0..dimensions[0]).map(|x| x as f32))
        .collect();
    let velocity_y = vec![0.0; config.element_count()];

    for (velocity, delta) in [(1.0, -0.25), (-1.0, 0.25)] {
        let velocity_x = vec![velocity; config.element_count()];
        let mut output = vec![0.0; config.element_count()];
        kernel
            .execute(&scalar, &velocity_x, &velocity_y, config, &mut output)
            .expect("expected value");

        let mut expected = scalar.clone();
        for y in 1..dimensions[1] - 1 {
            for x in 1..dimensions[0] - 1 {
                expected[y * dimensions[0] + x] += delta;
            }
        }
        assert_eq!(output, expected, "velocity_x={velocity}");
    }
}

#[test]
fn rejects_length_nonfinite_and_cfl_contract_violations() {
    let config = AdvectionConfig::new([3, 3, 1], [1.0; 3], 0.25).expect("expected value");
    let scalar = vec![1.0; config.element_count()];
    let velocity = vec![0.0; config.element_count()];
    let mut output = vec![0.0; config.element_count()];
    let Some(kernel) = kernel() else { return };

    let length_error = kernel
        .execute(&scalar[..8], &velocity, &velocity, config, &mut output)
        .expect_err("should reject insufficient scalar field size");
    assert!(matches!(
        length_error,
        Error::DimensionMismatch {
            expected: 9,
            actual: 8
        }
    ));

    let mut nonfinite_velocity = velocity.clone();
    nonfinite_velocity[4] = f32::NAN;
    let nonfinite_error = kernel
        .execute(&scalar, &nonfinite_velocity, &velocity, config, &mut output)
        .expect_err("should reject non-finite velocity values");
    assert!(matches!(nonfinite_error, Error::PhysicsViolation(_)));

    let unstable_velocity = vec![5.0; config.element_count()];
    let cfl_error = kernel
        .execute(&scalar, &unstable_velocity, &velocity, config, &mut output)
        .expect_err("should reject CFL-unstable velocity");
    assert!(matches!(cfl_error, Error::PhysicsViolation(_)));
}

#[test]
fn configuration_rejects_degenerate_grid_spacing_and_timestep() {
    assert!(matches!(
        AdvectionConfig::new([1, 3, 1], [1.0; 3], 0.1),
        Err(Error::InvalidConfiguration(_))
    ));
    assert!(matches!(
        AdvectionConfig::new([3, 3, 1], [0.0, 1.0, 1.0], 0.1),
        Err(Error::InvalidConfiguration(_))
    ));
    assert!(matches!(
        AdvectionConfig::new([3, 3, 1], [1.0; 3], -0.1),
        Err(Error::InvalidConfiguration(_))
    ));
}
