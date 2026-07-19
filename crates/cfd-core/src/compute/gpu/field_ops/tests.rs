use super::GpuFieldOps;
use crate::compute::gpu::GpuContext;
use crate::error::Error;
use aequitas::systems::si::{quantities::Length, units::Meter};
use std::sync::Arc;

fn operations() -> GpuFieldOps {
    let context = GpuContext::create().expect("GPU arithmetic tests require a WGPU provider");
    GpuFieldOps::new(Arc::new(context)).expect("field kernels must compile through Hephaestus")
}

fn meters(value: f32) -> Length<f32> {
    Length::from_unit::<Meter>(value)
}

#[test]
fn add_fields_handles_partial_workgroups() {
    let operations = operations();

    for len in [1, 63, 64, 65, 257] {
        let left: Vec<f32> = (0..len).map(|index| index as f32).collect();
        let right: Vec<f32> = (0..len).map(|index| 256.0 - index as f32).collect();
        let mut result = vec![0.0; len];

        operations
            .add_fields(&left, &right, &mut result)
            .expect("Hephaestus field addition must dispatch");

        assert_eq!(result, vec![256.0; len]);
    }
}

#[test]
fn multiply_field_handles_partial_workgroups() {
    let operations = operations();

    for len in [1, 63, 64, 65, 257] {
        let field: Vec<f32> = (0..len).map(|index| index as f32 * 0.5).collect();
        let expected: Vec<f32> = (0..len).map(|index| -(index as f32)).collect();
        let mut result = vec![0.0; len];

        operations
            .multiply_field(&field, -2.0, &mut result)
            .expect("Hephaestus scalar multiplication must dispatch");

        assert_eq!(result, expected);
    }
}

#[test]
fn arithmetic_accepts_empty_fields() {
    let operations = operations();
    let mut result = Vec::new();

    operations
        .add_fields(&[], &[], &mut result)
        .expect("empty addition is a valid no-op");
    operations
        .multiply_field(&[], 2.0, &mut result)
        .expect("empty multiplication is a valid no-op");

    assert!(result.is_empty());
}

#[test]
fn arithmetic_rejects_mismatched_lengths() {
    let operations = operations();
    let mut two_values = [0.0; 2];
    let mut one_value = [0.0; 1];

    let add_operand_error = operations
        .add_fields(&[1.0, 2.0], &[3.0], &mut two_values)
        .unwrap_err();
    assert!(matches!(
        add_operand_error,
        Error::DimensionMismatch {
            expected: 2,
            actual: 1
        }
    ));

    let add_output_error = operations
        .add_fields(&[1.0, 2.0], &[3.0, 4.0], &mut one_value)
        .unwrap_err();
    assert!(matches!(
        add_output_error,
        Error::DimensionMismatch {
            expected: 2,
            actual: 1
        }
    ));

    let multiply_output_error = operations
        .multiply_field(&[1.0, 2.0], 3.0, &mut one_value)
        .unwrap_err();
    assert!(matches!(
        multiply_output_error,
        Error::DimensionMismatch {
            expected: 2,
            actual: 1
        }
    ));
}

#[test]
fn laplacian_rejects_invalid_contracts() {
    let operations = operations();

    let mut output = [0.0; 4];
    let input_error = operations
        .laplacian_2d(&[1.0; 3], 2, 2, meters(1.0), meters(1.0), &mut output)
        .unwrap_err();
    assert!(matches!(
        input_error,
        Error::DimensionMismatch {
            expected: 4,
            actual: 3
        }
    ));

    let mut short_output = [0.0; 3];
    let output_error = operations
        .laplacian_2d(&[1.0; 4], 2, 2, meters(1.0), meters(1.0), &mut short_output)
        .unwrap_err();
    assert!(matches!(
        output_error,
        Error::DimensionMismatch {
            expected: 4,
            actual: 3
        }
    ));

    let spacing_error = operations
        .laplacian_2d(&[1.0; 4], 2, 2, meters(0.0), meters(1.0), &mut output)
        .unwrap_err();
    assert!(matches!(spacing_error, Error::InvalidConfiguration(_)));
}
