use super::{DiffusionConfig, GpuDiffusionKernel};
use crate::compute::gpu::GpuContext;
use crate::error::Error;
use std::sync::Arc;

fn kernel() -> GpuDiffusionKernel {
    let context = GpuContext::create().expect("diffusion tests require a WGPU provider");
    GpuDiffusionKernel::new(Arc::new(context))
        .expect("diffusion kernel must compile through Hephaestus")
}

#[test]
fn constant_field_is_exact_identity_across_partial_workgroups() {
    let config = DiffusionConfig::new([9, 5, 3], [1.0; 3], 0.125, 1.0).unwrap();
    let input = vec![7.25; config.element_count()];
    let mut output = vec![0.0; config.element_count()];

    kernel().execute(&input, config, &mut output).unwrap();

    assert_eq!(output, input);
}

#[test]
fn quadratic_field_has_exact_laplacian_and_copied_boundaries() {
    let dimensions = [5, 4, 3];
    let config = DiffusionConfig::new(dimensions, [1.0; 3], 0.125, 1.0).unwrap();
    let input: Vec<f32> = (0..dimensions[2])
        .flat_map(|z| {
            (0..dimensions[1])
                .flat_map(move |y| (0..dimensions[0]).map(move |x| (x * x + y * y + z * z) as f32))
        })
        .collect();
    let mut output = vec![0.0; config.element_count()];

    kernel().execute(&input, config, &mut output).unwrap();

    let mut expected = input.clone();
    for z in 1..dimensions[2] - 1 {
        for y in 1..dimensions[1] - 1 {
            for x in 1..dimensions[0] - 1 {
                let index = z * dimensions[0] * dimensions[1] + y * dimensions[0] + x;
                expected[index] += 0.75;
            }
        }
    }
    assert_eq!(output, expected);
}

#[test]
fn rejects_length_and_nonfinite_input() {
    let config = DiffusionConfig::new([3, 3, 3], [1.0; 3], 0.125, 1.0).unwrap();
    let input = vec![1.0; config.element_count()];
    let mut output = vec![0.0; config.element_count()];
    let kernel = kernel();

    let length_error = kernel
        .execute(&input[..26], config, &mut output)
        .unwrap_err();
    assert!(matches!(
        length_error,
        Error::DimensionMismatch {
            expected: 27,
            actual: 26
        }
    ));

    let mut nonfinite = input;
    nonfinite[13] = f32::NAN;
    let nonfinite_error = kernel.execute(&nonfinite, config, &mut output).unwrap_err();
    assert!(matches!(nonfinite_error, Error::PhysicsViolation(_)));
}

#[test]
fn configuration_rejects_invalid_and_unstable_parameters() {
    for result in [
        DiffusionConfig::new([2, 3, 3], [1.0; 3], 0.1, 1.0),
        DiffusionConfig::new([3, 3, 3], [0.0, 1.0, 1.0], 0.1, 1.0),
        DiffusionConfig::new([3, 3, 3], [1.0; 3], -0.1, 1.0),
        DiffusionConfig::new([3, 3, 3], [1.0; 3], 0.1, -1.0),
        DiffusionConfig::new([3, 3, 3], [1.0; 3], 0.25, 1.0),
    ] {
        assert!(matches!(result, Err(Error::InvalidConfiguration(_))));
    }
}
