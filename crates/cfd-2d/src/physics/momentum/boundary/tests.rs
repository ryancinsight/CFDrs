use super::*;
use crate::grid::StructuredGrid2D;
use cfd_core::physics::boundary::{BoundaryCondition, WallType};
use leto::geometry::Vector3;
use leto::Array1;
use std::collections::HashMap;

struct RecordingMatrix<T> {
    entries: Vec<(usize, usize, T)>,
}

impl<T> RecordingMatrix<T> {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }
}

impl<T: cfd_core::CfdScalar + Copy> MatrixUpdater<T> for RecordingMatrix<T> {
    fn add_entry(&mut self, row: usize, col: usize, val: T) -> cfd_core::error::Result<()> {
        self.entries.push((row, col, val));
        Ok(())
    }
}

#[test]
fn higher_order_west_boundary_skips_shared_corners() {
    let grid = StructuredGrid2D::new(4, 4, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)
        .expect("grid creation failed");
    let mut matrix = RecordingMatrix::new();
    let mut rhs = Array1::from_elem([grid.nx * grid.ny], 0.0);
    let boundaries = HashMap::from([(
        "west".to_string(),
        BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        },
    )]);

    apply_higher_order_wall_boundaries(
        &mut matrix,
        &mut rhs,
        MomentumComponent::V,
        &boundaries,
        &grid,
    )
    .expect("higher-order boundary application failed");

    let top_left = (grid.ny - 1) * grid.nx;
    let bottom_left = 0_usize;

    assert!(
        matrix.entries.iter().all(|(row, _, _)| *row != top_left),
        "west higher-order boundary must not write the top-left corner"
    );
    assert!(
        matrix.entries.iter().all(|(row, _, _)| *row != bottom_left),
        "west higher-order boundary must not write the bottom-left corner"
    );
    assert!(
        matrix.entries.iter().any(|(row, _, _)| *row == grid.nx),
        "west higher-order boundary must still write interior west-wall rows"
    );
}

#[test]
fn north_wall_owns_the_top_corners() {
    let grid = StructuredGrid2D::new(4, 4, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)
        .expect("grid creation failed");
    let mut matrix = RecordingMatrix::new();
    let mut rhs = Array1::from_elem([grid.nx * grid.ny], 0.0);
    let boundaries = HashMap::from([
        (
            "north".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::Moving {
                    velocity: Vector3::new(1.0, 0.0, 0.0),
                },
            },
        ),
        (
            "south".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        ),
        (
            "west".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        ),
        (
            "east".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        ),
    ]);

    apply_higher_order_wall_boundaries(
        &mut matrix,
        &mut rhs,
        MomentumComponent::U,
        &boundaries,
        &grid,
    )
    .expect("higher-order boundary application failed");
    apply_momentum_boundaries(
        &mut matrix,
        &mut rhs,
        MomentumComponent::U,
        &boundaries,
        &grid,
    )
    .expect("boundary application failed");

    let top_left = (grid.ny - 1) * grid.nx;
    let top_right = top_left + grid.nx - 1;
    let bottom_left = 0_usize;
    let bottom_right = grid.nx - 1;

    assert_eq!(rhs[top_left], 1.0);
    assert_eq!(rhs[top_right], 1.0);
    assert_eq!(rhs[bottom_left], 0.0);
    assert_eq!(rhs[bottom_right], 0.0);
}

/// The characteristic outlet writes the zero-gradient stencil, and nothing
/// else (CFDRS-OUTLET-FLAG-INERT).
///
/// This is the test whose absence let an inert flag survive: the variant once
/// carried an `extrapolate_velocity: bool` selecting between two branches that
/// applied the same three matrix entries, so no configuration of it could be
/// distinguished by behaviour. Asserting the entries rather than "it ran"
/// pins the one stencil the variant can express.
#[test]
fn characteristic_outlet_writes_the_zero_gradient_stencil() {
    let grid = StructuredGrid2D::new(4, 4, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)
        .expect("grid creation failed");
    let mut matrix = RecordingMatrix::new();
    let mut rhs = Array1::from_elem([grid.nx * grid.ny], 7.0_f64);
    let outlet = BoundaryCondition::characteristic_outlet(101_325.0_f64);

    apply_west_boundary(
        &mut matrix,
        &mut rhs,
        &outlet,
        MomentumComponent::U,
        &grid,
        grid.nx,
        grid.ny,
    )
    .expect("west outlet application");

    // Interior rows only: the horizontal walls own the corners, so a 4x4 grid
    // leaves rows j = 1 and j = 2 on the west edge.
    let expected: Vec<(usize, usize, f64)> = [1_usize, 2]
        .into_iter()
        .flat_map(|j| {
            let idx = j * grid.nx;
            [(idx, idx, 1.0_f64), (idx, idx + 1, -1.0_f64)]
        })
        .collect();
    assert_eq!(
        matrix.entries, expected,
        "the outlet must write u[boundary] - u[interior] = 0 and nothing else"
    );

    // The right-hand side is zeroed at those rows and untouched elsewhere, so
    // the pressure the variant carries does not leak into the momentum system.
    for j in 0..grid.ny {
        let idx = j * grid.nx;
        let expected_rhs = if (1..=2).contains(&j) { 0.0 } else { 7.0 };
        assert_eq!(rhs[idx], expected_rhs, "west row j = {j} right-hand side");
    }
}
