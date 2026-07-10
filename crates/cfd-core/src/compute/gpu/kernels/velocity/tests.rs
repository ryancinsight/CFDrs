use super::{GpuVelocityKernel, VelocityConfig};
use crate::compute::gpu::GpuContext;
use crate::error::Error;
use std::sync::Arc;

fn kernel() -> GpuVelocityKernel {
    let context = GpuContext::create().expect("velocity tests require a WGPU provider");
    GpuVelocityKernel::new(Arc::new(context))
        .expect("velocity kernels must compile through Hephaestus")
}

fn coordinates(dimensions: [usize; 3], field: impl Fn(usize, usize, usize) -> f32) -> Vec<f32> {
    (0..dimensions[2])
        .flat_map(|z| {
            let field = &field;
            (0..dimensions[1]).flat_map(move |y| (0..dimensions[0]).map(move |x| field(x, y, z)))
        })
        .collect()
}

fn interior_indices(dimensions: [usize; 3]) -> impl Iterator<Item = usize> {
    (1..dimensions[2] - 1).flat_map(move |z| {
        (1..dimensions[1] - 1).flat_map(move |y| {
            (1..dimensions[0] - 1)
                .map(move |x| z * dimensions[0] * dimensions[1] + y * dimensions[0] + x)
        })
    })
}

#[test]
fn correction_matches_linear_pressure_gradient_and_zeroes_boundaries() {
    let dimensions = [9, 5, 3];
    let config = VelocityConfig::new(dimensions, [1.0; 3], 0.5, 2.0).unwrap();
    let pressure = coordinates(dimensions, |x, y, z| {
        (2 * x) as f32 - (4 * y) as f32 + (6 * z) as f32
    });
    let velocity_x = vec![10.0; config.element_count()];
    let velocity_y = vec![20.0; config.element_count()];
    let velocity_z = vec![30.0; config.element_count()];
    let mut output_x = vec![-1.0; config.element_count()];
    let mut output_y = vec![-1.0; config.element_count()];
    let mut output_z = vec![-1.0; config.element_count()];

    kernel()
        .correct(
            [&velocity_x, &velocity_y, &velocity_z],
            &pressure,
            config,
            [&mut output_x, &mut output_y, &mut output_z],
        )
        .unwrap();

    let mut expected_x = vec![0.0; config.element_count()];
    let mut expected_y = vec![0.0; config.element_count()];
    let mut expected_z = vec![0.0; config.element_count()];
    for index in interior_indices(dimensions) {
        expected_x[index] = 9.5;
        expected_y[index] = 21.0;
        expected_z[index] = 28.5;
    }
    assert_eq!(output_x, expected_x);
    assert_eq!(output_y, expected_y);
    assert_eq!(output_z, expected_z);
}

#[test]
fn divergence_source_matches_linear_velocity_field() {
    let dimensions = [5, 4, 3];
    let config = VelocityConfig::new(dimensions, [1.0; 3], 0.5, 2.0).unwrap();
    let velocity_x = coordinates(dimensions, |x, _, _| (2 * x) as f32);
    let velocity_y = coordinates(dimensions, |_, y, _| -(4 * y as i32) as f32);
    let velocity_z = coordinates(dimensions, |_, _, z| (6 * z) as f32);
    let mut output = vec![-1.0; config.element_count()];

    kernel()
        .divergence_source(&velocity_x, &velocity_y, &velocity_z, config, &mut output)
        .unwrap();

    let mut expected = vec![0.0; config.element_count()];
    for index in interior_indices(dimensions) {
        expected[index] = 16.0;
    }
    assert_eq!(output, expected);
}

#[test]
fn rejects_length_and_nonfinite_fields() {
    let config = VelocityConfig::new([3, 3, 3], [1.0; 3], 0.5, 2.0).unwrap();
    let field = vec![1.0; config.element_count()];
    let mut output = vec![0.0; config.element_count()];
    let kernel = kernel();

    let length_error = kernel
        .divergence_source(&field[..26], &field, &field, config, &mut output)
        .unwrap_err();
    assert!(matches!(
        length_error,
        Error::DimensionMismatch {
            expected: 27,
            actual: 26
        }
    ));

    let mut nonfinite = field.clone();
    nonfinite[13] = f32::INFINITY;
    let nonfinite_error = kernel
        .divergence_source(&field, &nonfinite, &field, config, &mut output)
        .unwrap_err();
    assert!(matches!(nonfinite_error, Error::PhysicsViolation(_)));
}

#[test]
fn configuration_rejects_degenerate_and_nonphysical_parameters() {
    for result in [
        VelocityConfig::new([2, 3, 3], [1.0; 3], 0.5, 2.0),
        VelocityConfig::new([3, 3, 3], [0.0, 1.0, 1.0], 0.5, 2.0),
        VelocityConfig::new([3, 3, 3], [1.0; 3], 0.0, 2.0),
        VelocityConfig::new([3, 3, 3], [1.0; 3], 0.5, 0.0),
    ] {
        assert!(matches!(result, Err(Error::InvalidConfiguration(_))));
    }
}
