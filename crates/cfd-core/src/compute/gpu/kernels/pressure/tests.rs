use super::{GpuPressureKernel, PressureConfig};
use crate::compute::gpu::GpuContext;
use crate::error::Error;
use std::sync::Arc;

fn kernel() -> GpuPressureKernel {
    let context = GpuContext::create().expect("pressure tests require a WGPU provider");
    GpuPressureKernel::new(Arc::new(context))
        .expect("pressure kernels must compile through Hephaestus")
}

fn coordinates(dimensions: [usize; 3], field: impl Fn(usize, usize, usize) -> f32) -> Vec<f32> {
    (0..dimensions[2])
        .flat_map(|z| {
            let field = &field;
            (0..dimensions[1]).flat_map(move |y| (0..dimensions[0]).map(move |x| field(x, y, z)))
        })
        .collect()
}

fn index(dimensions: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    z * dimensions[0] * dimensions[1] + y * dimensions[0] + x
}

#[test]
fn iteration_preserves_quadratic_solution_and_clamps_neumann_boundaries() {
    let dimensions = [9, 5, 3];
    let config = PressureConfig::new(dimensions, [1.0; 3], 0.5).unwrap();
    let pressure = coordinates(dimensions, |x, y, z| (x * x + y * y + z * z) as f32);
    let source = vec![6.0; config.element_count()];
    let mut output = vec![-1.0; config.element_count()];

    kernel()
        .iterate(&pressure, &source, config, &mut output)
        .unwrap();

    let expected = coordinates(dimensions, |x, y, z| {
        let x = x.clamp(1, dimensions[0] - 2);
        let y = y.clamp(1, dimensions[1] - 2);
        let z = z.clamp(1, dimensions[2] - 2);
        (x * x + y * y + z * z) as f32
    });
    assert_eq!(output, expected);
}

#[test]
fn iteration_applies_weighted_source_term_exactly() {
    let dimensions = [3, 3, 3];
    let config = PressureConfig::new(dimensions, [1.0; 3], 0.5).unwrap();
    let pressure = vec![0.0; config.element_count()];
    let source = vec![-12.0; config.element_count()];
    let mut output = vec![-1.0; config.element_count()];

    kernel()
        .iterate(&pressure, &source, config, &mut output)
        .unwrap();

    let mut expected = vec![0.0; config.element_count()];
    expected[index(dimensions, 1, 1, 1)] = 1.0;
    assert_eq!(output, expected);
}

#[test]
fn residual_matches_quadratic_laplacian_and_zeroes_boundaries() {
    let dimensions = [5, 4, 3];
    let config = PressureConfig::new(dimensions, [1.0; 3], 1.0).unwrap();
    let pressure = coordinates(dimensions, |x, y, z| (x * x + y * y + z * z) as f32);
    let source = vec![4.0; config.element_count()];
    let mut output = vec![-1.0; config.element_count()];

    kernel()
        .residual(&pressure, &source, config, &mut output)
        .unwrap();

    let mut expected = vec![0.0; config.element_count()];
    for z in 1..dimensions[2] - 1 {
        for y in 1..dimensions[1] - 1 {
            for x in 1..dimensions[0] - 1 {
                expected[index(dimensions, x, y, z)] = 2.0;
            }
        }
    }
    assert_eq!(output, expected);
}

#[test]
fn rejects_length_and_nonfinite_fields() {
    let config = PressureConfig::new([3, 3, 3], [1.0; 3], 1.0).unwrap();
    let field = vec![1.0; config.element_count()];
    let mut output = vec![0.0; config.element_count()];
    let kernel = kernel();

    let length_error = kernel
        .residual(&field[..26], &field, config, &mut output)
        .unwrap_err();
    assert!(matches!(
        length_error,
        Error::DimensionMismatch {
            expected: 27,
            actual: 26
        }
    ));

    let mut nonfinite = field.clone();
    nonfinite[13] = f32::NAN;
    let nonfinite_error = kernel
        .residual(&field, &nonfinite, config, &mut output)
        .unwrap_err();
    assert!(matches!(nonfinite_error, Error::PhysicsViolation(_)));
}

#[test]
fn configuration_rejects_degenerate_and_invalid_parameters() {
    for result in [
        PressureConfig::new([2, 3, 3], [1.0; 3], 1.0),
        PressureConfig::new([3, 3, 3], [0.0, 1.0, 1.0], 1.0),
        PressureConfig::new([3, 3, 3], [1.0; 3], 0.0),
        PressureConfig::new([3, 3, 3], [1.0; 3], 1.1),
    ] {
        assert!(matches!(result, Err(Error::InvalidConfiguration(_))));
    }
}
