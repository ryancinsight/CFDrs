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
            tracing::debug!(error = ?error, "GPU turbulence test unavailable");
            None
        }
    }
}

#[test]
fn smagorinsky_matches_linear_strain_and_zeroes_boundaries() {
    let Some(compute) = compute() else { return };
    let dimensions = [9, 5];
    let grid = TurbulenceGrid::new(dimensions, [1.0; 2]).expect("invariant: valid turbulence grid");
    let velocity_x = coordinates(dimensions, |x, _| x as f32);
    let velocity_y = coordinates(dimensions, |_, y| y as f32);
    let mut output = vec![-1.0; grid.element_count()];

    compute
        .compute_smagorinsky_sgs(&velocity_x, &velocity_y, grid, 0.5, &mut output)
        .expect("invariant: valid Smagorinsky computation");

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
    let grid = TurbulenceGrid::new([5, 3], [0.25, 1.0]).expect("invariant: valid turbulence grid");
    let mut output = vec![0.0; grid.element_count()];

    compute
        .compute_des_length_scale(grid, 0.5, &mut output)
        .expect("invariant: valid DES computation");

    assert_eq!(output, vec![0.25; grid.element_count()]);
}

#[test]
fn wall_distance_matches_rectangular_geometry() {
    let Some(compute) = compute() else { return };
    let dimensions = [5, 5];
    let grid =
        TurbulenceGrid::new(dimensions, [0.5, 0.25]).expect("invariant: valid turbulence grid");
    let mut output = vec![-1.0; grid.element_count()];

    compute
        .compute_wall_distance(grid, &mut output)
        .expect("invariant: valid wall-distance computation");

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
    let grid = TurbulenceGrid::new([3, 3], [1.0; 2]).expect("invariant: valid turbulence grid");
    let field = vec![0.0; grid.element_count()];
    let mut output = vec![0.0; grid.element_count()];
    let Some(compute) = compute() else { return };
    let length_error = compute
        .compute_smagorinsky_sgs(&field[..8], &field, grid, 0.1, &mut output)
        .expect_err("invariant: short turbulence field is rejected");
    assert!(matches!(length_error, Error::DimensionMismatch { .. }));
    let constant_error = compute
        .compute_des_length_scale(grid, -0.1, &mut output)
        .expect_err("invariant: negative DES scale is rejected");
    assert!(matches!(constant_error, Error::InvalidConfiguration(_)));
    let mut nonfinite = field.clone();
    nonfinite[4] = f32::NAN;
    let value_error = compute
        .compute_smagorinsky_sgs(&field, &nonfinite, grid, 0.1, &mut output)
        .expect_err("invariant: nonfinite turbulence field is rejected");
    assert!(matches!(value_error, Error::PhysicsViolation(_)));
}
