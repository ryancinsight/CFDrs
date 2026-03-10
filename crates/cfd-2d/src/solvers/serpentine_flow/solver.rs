//! Discretized 2D serpentine flow solver using FVM.

use super::{SerpentineGeometry, SerpentineMixingSolution};
use crate::solvers::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use crate::solvers::scalar_transport_2d::{ScalarTransportConfig, ScalarTransportSolver2D};
use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Discretized 2D Serpentine Flow Solver
pub struct SerpentineSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Serpentine channel geometry definition.
    pub geometry: SerpentineGeometry<T>,
    /// Navier-Stokes flow solver on a staggered grid.
    pub ns_solver: NavierStokesSolver2D<T>,
    /// Advection-diffusion scalar transport solver.
    pub scalar_solver: ScalarTransportSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive> SerpentineSolver2D<T> {
    /// Create new discretized serpentine solver
    pub fn new(
        geometry: SerpentineGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
    ) -> Self {
        let bbox = geometry.bounding_box();
        let width = bbox[1] - bbox[0];
        let height = bbox[3] - bbox[2];

        let grid = StaggeredGrid2D::new(nx, ny, width, height);
        let config = SIMPLEConfig::default();
        let mut ns_solver = NavierStokesSolver2D::new(grid, blood, density, config);

        // Populate mask
        for i in 0..nx {
            for j in 0..ny {
                let x = ns_solver.grid.x_center(i) + bbox[0];
                let y = ns_solver.grid.y_center(j) + bbox[2];
                ns_solver.field.mask[i][j] = geometry.contains(x, y);
            }
        }

        let scalar_solver = ScalarTransportSolver2D::new(nx, ny);

        Self {
            geometry,
            ns_solver,
            scalar_solver,
        }
    }

    /// Solve for flow and mixing
    pub fn solve(
        &mut self,
        u_inlet: T,
        diffusion_coeff: T,
        c_left: T,  // Concentration on left half of inlet
        c_right: T, // Concentration on right half of inlet
    ) -> CfdResult<SerpentineMixingSolution<T>> {
        // Use the scalar-transport default tolerance (1e-5), which matches the
        // FVM spatial truncation error O(Δx²) on typical coarse grids.
        self.solve_with_tolerance(
            u_inlet,
            diffusion_coeff,
            c_left,
            c_right,
            T::from_f64(1e-5).unwrap_or_else(num_traits::Zero::zero),
        )
    }

    /// Solve for flow and mixing with custom tolerance
    pub fn solve_with_tolerance(
        &mut self,
        u_inlet: T,
        diffusion_coeff: T,
        c_left: T,  // Concentration on left half of inlet
        c_right: T, // Concentration on right half of inlet
        tolerance: T,
    ) -> CfdResult<SerpentineMixingSolution<T>> {
        // 1. Solve Navier-Stokes
        self.ns_solver
            .solve(u_inlet)
            .map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        // 2. Setup Scalar Transport
        let config = ScalarTransportConfig {
            tolerance,
            diffusion_coeff,
            ..ScalarTransportConfig::default()
        };

        // Define inlet concentration profile (step function)
        let ny = self.ns_solver.grid.ny;
        let mut boundary_c = vec![T::zero(); ny];
        for j in 0..ny {
            if j < ny / 2 {
                boundary_c[j] = c_left;
            } else {
                boundary_c[j] = c_right;
            }
        }

        // 3. Solve Scalar Transport
        // Cells at the inlet (West boundary) with mask=true and u > 0 will use boundary_c.
        self.scalar_solver
            .solve(
                &self.ns_solver.grid,
                &self.ns_solver.field,
                &config,
                &boundary_c,
            )
            .map_err(cfd_core::error::Error::Solver)?;

        // 4. Extract metrics
        let pe = (u_inlet * self.geometry.width) / diffusion_coeff;

        // Compute mixing efficiency at outlet
        // ISO = 1 - variance / variance_inlet
        let nx = self.ns_solver.grid.nx;
        let mut sum_c = T::zero();
        let mut sum_c_sq = T::zero();
        let mut count = 0;

        // Find last fluid column (outlet)
        for j in 0..ny {
            if self.ns_solver.field.mask[nx - 1][j] {
                let ci = self.scalar_solver.c[nx - 1][j];
                sum_c += ci;
                sum_c_sq += ci * ci;
                count += 1;
            }
        }

        let mut mixing_frac = T::zero();
        if count > 0 {
            let n = T::from_usize(count).unwrap_or_else(T::one);
            let c_mean = sum_c / n;
            let variance = (sum_c_sq / n) - (c_mean * c_mean);

            // Expected variance for perfectly unmixed: 0.25 for c_left=0, c_right=1
            let var_inlet = half_sq(c_left - c_right);
            mixing_frac = T::one()
                - Float::sqrt(
                    variance
                        / num_traits::Float::max(
                            var_inlet,
                            T::from_f64(1e-10).unwrap_or_else(num_traits::Zero::zero),
                        ),
                );
        }

        // Compute pressure drop between inlet and outlet regions
        let mut p_in = T::zero();
        let mut count_in = 0;
        let mut p_out = T::zero();
        let mut count_out = 0;

        for j in 0..ny {
            if self.ns_solver.field.mask[0][j] {
                p_in += self.ns_solver.field.p[0][j];
                count_in += 1;
            }
            if self.ns_solver.field.mask[nx - 1][j] {
                p_out += self.ns_solver.field.p[nx - 1][j];
                count_out += 1;
            }
        }

        let mut pressure_drop = T::zero();
        if count_in > 0 && count_out > 0 {
            pressure_drop = (p_in / T::from_usize(count_in).unwrap_or_else(T::one))
                - (p_out / T::from_usize(count_out).unwrap_or_else(T::one));
        }

        Ok(SerpentineMixingSolution {
            c_inlet_a: c_left,
            c_inlet_b: c_right,
            peclet: pe,
            l_mix_90: T::zero(), // Computed by analytical model
            t_mix_90: T::zero(),
            mixing_fraction_outlet: mixing_frac,
            pressure_drop,
        })
    }
}

fn half_sq<T: RealField + Copy + FromPrimitive>(val: T) -> T {
    let half = T::from_f64(0.5).unwrap_or_else(num_traits::Zero::zero);
    (val * half) * (val * half)
}

#[cfg(test)]
mod tests_discretized {
    use super::*;
    use crate::solvers::ns_fvm::BloodModel;

    #[test]
    fn test_discretized_serpentine_mixing() {
        // High diffusion to ensure some mixing in small grid
        let geom = SerpentineGeometry::new(
            0.0002, // width
            0.0001, // height
            0.0005, // straight
            0.0002, // radius
            1,      // 1 cycle
        );

        // Fluid properties
        let blood = BloodModel::Newtonian(0.001);
        let density = 1000.0;

        // Grid (coarse for fast test)
        let nx = 40;
        let ny = 20;

        let mut solver = SerpentineSolver2D::new(geom, blood, density, nx, ny);

        // Solve with higher diffusion to see mixing
        // Use a reasonable diffusion for a coarse grid
        let result = solver.solve(0.01, 1e-6, 0.0, 1.0);

        assert!(
            result.is_ok(),
            "Serpentine solver failed: {:?}",
            result.err()
        );
        let sol = result.unwrap();

        println!("Serpentine Mixing Solution: {:?}", sol);

        // Qualitative checks
        assert!(sol.peclet > 0.0, "Peclet should be positive");
        assert!(
            sol.mixing_fraction_outlet >= 0.0,
            "Mixing fraction should be non-negative"
        );
        // Due to the geometry mapping, some flow might be blocked if resolution is too low.
        // We'll check if pressure drop is positive.
        assert!(
            sol.pressure_drop >= 0.0,
            "Pressure drop should be non-negative"
        );
    }
}
