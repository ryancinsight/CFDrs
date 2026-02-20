//! [`ReynoldsStressModel`] — main model struct, constructor, wall-distance utilities,
//! and tensor initialisation.

use super::PressureStrainModel;
use super::tensor::ReynoldsStressTensor;
use super::super::constants::C_MU;
use nalgebra::{DMatrix, RealField};
use num_traits::{FromPrimitive, ToPrimitive};

/// Reynolds Stress Transport Model — model parameters and grid metadata.
#[derive(Debug, Clone)]
pub struct ReynoldsStressModel<T: RealField + Copy + FromPrimitive + ToPrimitive> {
    pub(super) nx: usize,
    pub(super) ny: usize,
    pub(super) wall_distance: Option<Vec<T>>,
    pub(super) c_mu: T,
    pub c1: T,
    pub(super) c2: T,
    pub(super) c1_star: T,
    pub(super) c2_star: T,
    pub(super) c3: T,
    pub(super) c3_star: T,
    pub(super) c4: T,
    pub(super) c5: T,
    pub pressure_strain_model: PressureStrainModel,
    pub(super) wall_reflection: bool,
    pub(super) curvature_correction: bool,
}

fn constant<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("RSTM constant must be representable")
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> ReynoldsStressModel<T> {
    /// Construct a new RSTM with literature-calibrated default constants.
    ///
    /// # Panics
    /// Panics if `nx == 0` or `ny == 0`.
    pub fn new(nx: usize, ny: usize) -> Self {
        assert!(nx > 0 && ny > 0, "Grid dimensions must be positive: nx={nx}, ny={ny}");
        Self {
            nx,
            ny,
            wall_distance: None,
            c_mu: constant(C_MU),
            c1: constant(1.8),
            c2: constant(0.6),
            c1_star: constant(1.7),
            c2_star: constant(-1.05),
            c3: constant(0.8),
            c3_star: constant(1.3),
            c4: constant(1.25),
            c5: constant(0.40),
            pressure_strain_model: PressureStrainModel::Quadratic,
            wall_reflection: true,
            curvature_correction: true,
        }
    }

    /// Inject a pre-computed wall-distance vector (nx × ny, row-major).
    pub fn set_wall_distance(&mut self, wall_distance: Vec<T>) {
        assert_eq!(
            wall_distance.len(), self.nx * self.ny,
            "Wall distance vector must have {} elements", self.nx * self.ny
        );
        self.wall_distance = Some(wall_distance);
    }

    /// Compute wall-distance field via iterative distance transform.
    pub fn compute_wall_distance_field(&self, wall_mask: &[bool], dx: T, dy: T) -> Vec<T> {
        let n = self.nx * self.ny;
        let max_val = T::max_value().expect("max_value required");
        let mut distance = vec![max_val; n];
        for i in 0..n { if wall_mask[i] { distance[i] = T::zero(); } }

        let mut changed = true;
        while changed {
            changed = false;
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let idx = y * self.nx + x;
                    if wall_mask[idx] { continue; }
                    let mut min_d = distance[idx];
                    let candidates = [
                        x.checked_sub(1).map(|xn| (y * self.nx + xn, dx)),
                        (x + 1 < self.nx).then_some((y * self.nx + x + 1, dx)),
                        y.checked_sub(1).map(|yn| (yn * self.nx + x, dy)),
                        (y + 1 < self.ny).then_some(((y + 1) * self.nx + x, dy)),
                    ];
                    for c in candidates.into_iter().flatten() {
                        let candidate = distance[c.0] + c.1;
                        if candidate < min_d { min_d = candidate; changed = true; }
                    }
                    distance[idx] = min_d;
                }
            }
        }
        distance
    }

    /// Initialise a [`ReynoldsStressTensor`] with isotropic turbulence.
    pub fn initialize_reynolds_stresses(&self, initial_k: T, initial_epsilon: T) -> ReynoldsStressTensor<T> {
        assert!(initial_k > T::zero(), "Initial TKE must be positive");
        assert!(initial_epsilon > T::zero(), "Initial dissipation must be positive");

        let iso = initial_k * constant::<T>(2.0 / 3.0);
        let mut xx = DMatrix::zeros(self.nx, self.ny);
        let xy = DMatrix::zeros(self.nx, self.ny);
        let mut yy = DMatrix::zeros(self.nx, self.ny);
        let mut k   = DMatrix::zeros(self.nx, self.ny);
        let mut eps = DMatrix::zeros(self.nx, self.ny);

        for i in 0..self.nx { for j in 0..self.ny {
            xx[(i, j)] = iso; yy[(i, j)] = iso;
            k[(i, j)] = initial_k; eps[(i, j)] = initial_epsilon;
        } }

        ReynoldsStressTensor { xx, xy, yy, k, epsilon: eps, epsilon_xx: None, epsilon_xy: None, epsilon_yy: None }
    }

    /// Activate anisotropic dissipation tensor components.
    pub fn enable_dissipation_tensor(&self, tensor: &mut ReynoldsStressTensor<T>, initial_epsilon: T) {
        let iso = initial_epsilon * constant::<T>(2.0 / 3.0);
        let mut eps_xx = DMatrix::zeros(self.nx, self.ny);
        let eps_xy = DMatrix::zeros(self.nx, self.ny);
        let mut eps_yy = DMatrix::zeros(self.nx, self.ny);
        for i in 0..self.nx { for j in 0..self.ny {
            eps_xx[(i, j)] = iso; eps_yy[(i, j)] = iso;
        } }
        tensor.epsilon_xx = Some(eps_xx);
        tensor.epsilon_xy = Some(eps_xy);
        tensor.epsilon_yy = Some(eps_yy);
    }
}
