#![allow(clippy::float_cmp)]
#![allow(clippy::print_stdout)]
#![allow(clippy::print_stderr)]
use super::*;
use crate::error::Error;

fn coordinates(dimensions: [usize; 2], field: impl Fn(usize, usize) -> f32) -> Vec<f32> {
    (0..dimensions[1])
        .flat_map(|y| {
            (0..dimensions[0]).map({
                let field = &field;
                move |x| field(x, y)
            })
        })
        .collect()
}

fn compute() -> Option<GpuTurbulenceCompute> {
    match GpuTurbulenceCompute::new() {
        Ok(compute) => Some(compute),
        Err(error) => {
            eprintln!("Skipping GPU turbulence test: {error:?}");
            None
        }
    }
}

#[test]
fn smagorinsky_matches_linear_strain_and_zeroes_boundaries() {
    let Some(compute) = compute() else { return };
    let dimensions = [9, 5];
    let grid = TurbulenceGrid::new(dimensions, [1.0; 2]).expect("expected value");
    let velocity_x = coordinates(dimensions, |x, _| x as f32);
    let velocity_y = coordinates(dimensions, |_, y| y as f32);
    let mut output = vec![-1.0; grid.element_count()];

    compute
        .compute_smagorinsky_sgs(&velocity_x, &velocity_y, grid, 0.5, &mut output)
        .expect("expected value");

    let expected = coordinates(dimensions, |x, y| {
        if x == 0 || x == dimensions[0] - 1 || y == 0 || y == dimensions[1] - 1 {
            0.0
        } else {
            0.5
        }
    });
    assert_eq!(output, expected);
}

#[test]
fn des_grid_scale_is_input_independent_constant() {
    let Some(compute) = compute() else { return };
    let grid = TurbulenceGrid::new([5, 3], [0.25, 1.0]).expect("expected value");
    let mut output = vec![0.0; grid.element_count()];

    compute
        .compute_des_length_scale(grid, 0.5, &mut output)
        .expect("expected value");

    assert_eq!(output, vec![0.25; grid.element_count()]);
}

#[test]
fn wall_distance_matches_rectangular_geometry() {
    let Some(compute) = compute() else { return };
    let dimensions = [5, 5];
    let grid = TurbulenceGrid::new(dimensions, [0.5, 0.25]).expect("expected value");
    let mut output = vec![-1.0; grid.element_count()];

    compute
        .compute_wall_distance(grid, &mut output)
        .expect("expected value");

    let expected = coordinates(dimensions, |x, y| {
        let distance_x = x.min(dimensions[0] - 1 - x) as f32 * 0.5;
        let distance_y = y.min(dimensions[1] - 1 - y) as f32 * 0.25;
        distance_x.min(distance_y)
    });
    assert_eq!(output, expected);
}

#[test]
fn rejects_invalid_grid_lengths_constants_and_values() {
    assert!(matches!(
        TurbulenceGrid::new([2, 3], [1.0; 2]),
        Err(Error::InvalidConfiguration(_))
    ));
    let grid = TurbulenceGrid::new([3, 3], [1.0; 2]).expect("expected value");
    let field = vec![0.0; grid.element_count()];
    let mut output = vec![0.0; grid.element_count()];
    let Some(compute) = compute() else { return };
    let length_error = compute
        .compute_smagorinsky_sgs(&field[..8], &field, grid, 0.1, &mut output)
        .expect_err("should reject insufficient field size");
    assert!(matches!(length_error, Error::DimensionMismatch { .. }));
    let constant_error = compute
        .compute_des_length_scale(grid, -0.1, &mut output)
        .expect_err("should reject negative constant");
    assert!(matches!(constant_error, Error::InvalidConfiguration(_)));
    let mut nonfinite = field.clone();
    nonfinite[4] = f32::NAN;
    let value_error = compute
        .compute_smagorinsky_sgs(&field, &nonfinite, grid, 0.1, &mut output)
        .expect_err("should reject non-finite values");
    assert!(matches!(value_error, Error::PhysicsViolation(_)));
}
