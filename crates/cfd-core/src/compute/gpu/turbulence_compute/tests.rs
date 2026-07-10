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

#[test]
fn smagorinsky_matches_linear_strain_and_zeroes_boundaries() {
    let dimensions = [9, 5];
    let grid = TurbulenceGrid::new(dimensions, [1.0; 2]).unwrap();
    let velocity_x = coordinates(dimensions, |x, _| x as f32);
    let velocity_y = coordinates(dimensions, |_, y| y as f32);
    let mut output = vec![-1.0; grid.element_count()];

    GpuTurbulenceCompute::new()
        .unwrap()
        .compute_smagorinsky_sgs(&velocity_x, &velocity_y, grid, 0.5, &mut output)
        .unwrap();

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
    let grid = TurbulenceGrid::new([5, 3], [0.25, 1.0]).unwrap();
    let mut output = vec![0.0; grid.element_count()];

    GpuTurbulenceCompute::new()
        .unwrap()
        .compute_des_length_scale(grid, 0.5, &mut output)
        .unwrap();

    assert_eq!(output, vec![0.25; grid.element_count()]);
}

#[test]
fn wall_distance_matches_rectangular_geometry() {
    let dimensions = [5, 5];
    let grid = TurbulenceGrid::new(dimensions, [0.5, 0.25]).unwrap();
    let mut output = vec![-1.0; grid.element_count()];

    GpuTurbulenceCompute::new()
        .unwrap()
        .compute_wall_distance(grid, &mut output)
        .unwrap();

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
    let grid = TurbulenceGrid::new([3, 3], [1.0; 2]).unwrap();
    let field = vec![0.0; grid.element_count()];
    let mut output = vec![0.0; grid.element_count()];
    let compute = GpuTurbulenceCompute::new().unwrap();
    let length_error = compute
        .compute_smagorinsky_sgs(&field[..8], &field, grid, 0.1, &mut output)
        .unwrap_err();
    assert!(matches!(length_error, Error::DimensionMismatch { .. }));
    let constant_error = compute
        .compute_des_length_scale(grid, -0.1, &mut output)
        .unwrap_err();
    assert!(matches!(constant_error, Error::InvalidConfiguration(_)));
    let mut nonfinite = field.clone();
    nonfinite[4] = f32::NAN;
    let value_error = compute
        .compute_smagorinsky_sgs(&field, &nonfinite, grid, 0.1, &mut output)
        .unwrap_err();
    assert!(matches!(value_error, Error::PhysicsViolation(_)));
}
