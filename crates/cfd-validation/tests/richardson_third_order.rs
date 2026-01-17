use cfd_validation::geometry::RectangularDomain;
use cfd_validation::manufactured::richardson::MmsRichardsonStudy;
use cfd_validation::manufactured::ManufacturedDiffusion;

#[test]
fn richardson_estimates_third_order_uniform_ratio() -> anyhow::Result<()> {
    let mms = ManufacturedDiffusion::<f64>::new(1.0);
    let geometry = RectangularDomain::new(0.0, 1.0, 0.0, 1.0);

    let study = MmsRichardsonStudy::with_geometric_refinement(
        Box::new(mms),
        Box::new(geometry),
        4,   // number of grid levels (coarse to fine)
        1.0, // base grid size (unused by this test path)
        0.0, // evaluation time
    )?;

    // Use uniform refinement ratios: coarse→fine as 2, 4, 8, 16
    let grid_sizes = vec![2usize, 4, 8, 16];

    // Manufactured convergence: f(h) = h^3 (third-order)
    let solution_computer = |size: usize| {
        let h = 1.0 / size as f64;
        h.powi(3)
    };

    let result = study.richardson_extrapolation_error(&grid_sizes, solution_computer);

    // Verify the estimated order is ~3.0
    let p = result.estimated_order;
    assert!(p > 2.8 && p < 3.2, "Estimated order not ~3.0, got {}", p);

    // Sanity: extrapolated solution must be non-negative and finite
    let phi = result.extrapolated_solution;
    assert!(phi >= 0.0 && phi.is_finite());

    Ok(())
}

#[test]
fn richardson_estimates_third_order_nonuniform_ratio() -> anyhow::Result<()> {
    // Non-uniform refinement ratios to exercise bisection path: 12→6→3→2 (r21=2, r32=2, r43=1.5)
    let mms = ManufacturedDiffusion::<f64>::new(1.0);
    let geometry = RectangularDomain::new(0.0, 1.0, 0.0, 1.0);

    let study = MmsRichardsonStudy::with_geometric_refinement(
        Box::new(mms),
        Box::new(geometry),
        4,
        1.0,
        0.0,
    )?;

    // Use mildly nonuniform sizes to avoid perfect uniformity: 2, 3, 6, 12
    let grid_sizes = vec![2usize, 3, 6, 12];

    let solution_computer = |size: usize| {
        let h = 1.0 / size as f64;
        h.powi(3)
    };

    let result = study.richardson_extrapolation_error(&grid_sizes, solution_computer);

    let p = result.estimated_order;
    assert!(
        p > 2.6 && p < 3.4,
        "Estimated order not near 3.0 with nonuniform ratios, got {}",
        p
    );

    Ok(())
}
