//! Macroscopic quantity extraction for D2Q9 LBM.
//!
//! Computes density ρ, velocity **u**, and pressure *p* from the flat
//! distribution buffer `f[j * nx * 9 + i * 9 + q]`.
//!
//! # Theorem — Moment Consistency
//!
//! **Statement**: For any distribution $\{f_q\}$ the following zeroth and
//! first moments are exact (no approximation):
//!
//! $$\rho = \sum_{q=0}^{8} f_q, \qquad \rho \mathbf{u} = \sum_{q=0}^{8} f_q \mathbf{e}_q$$
//!
//! **Proof**:
//!
//! 1. The mapping $(f_q) \mapsto (\rho, \mathbf{j})$ is a linear projection whose
//!    correctness follows directly from the definitions; no truncation is performed.
//! 2. At equilibrium $f_q = f_q^{eq}$, the Chapman-Enskog analysis shows
//!    $\rho^{eq} = \rho$ and $\mathbf{j}^{eq} = \rho\mathbf{u}$ — the moments are
//!    consistent with the macroscopic fields by construction. □
//!
//! # Theorem — LBM Pressure Relation
//!
//! **Statement**: The LBM pressure is $p = c_s^2 \rho$ with $c_s^2 = 1/3$
//! (in lattice units).
//!
//! **Proof sketch**: From the Chapman-Enskog expansion at $O(\epsilon^0)$,
//! the leading-order equation of state is the ideal gas law $p = \rho c_s^2$
//! where $c_s = 1/\sqrt{3}$ is the lattice sound speed. This is exact
//! within the weakly-compressible limit $Ma \ll 1$. □

use crate::solvers::lbm::lattice::D2Q9;
use crate::solvers::lbm::streaming::f_idx;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Container for macroscopic fields stored contiguously.
///
/// # Layout
/// - `density[j * nx + i]`
/// - `velocity[(j * nx + i) * 2 + d]` (d=0 → x, d=1 → y)
/// - `pressure[j * nx + i]` (optional, only if `with_pressure()` called)
/// - `nuclei_fraction[j * nx + i]` (optional, only if `with_nuclei()` called)
#[derive(Debug, Clone)]
pub struct MacroscopicQuantities<T: RealField + Copy> {
    /// Flat density field; index `j * nx + i`.
    pub density: Vec<T>,
    /// Flat 2-component velocity field; index `(j * nx + i) * 2 + component`.
    pub velocity: Vec<T>,
    /// Optional flat pressure field; index `j * nx + i`.
    pub pressure: Option<Vec<T>>,
    /// Optional flat nuclei volume fraction field; index `j * nx + i`.
    pub nuclei_fraction: Option<Vec<T>>,
    /// Grid columns
    pub nx: usize,
    /// Grid rows
    pub ny: usize,
}

impl<T: RealField + Copy + FromPrimitive> MacroscopicQuantities<T> {
    /// Construct, initialising ρ = 1, **u** = 0 everywhere.
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            density: vec![T::one(); n],
            velocity: vec![T::zero(); n * 2],
            pressure: None,
            nuclei_fraction: None,
            nx,
            ny,
        }
    }

    /// Enable pressure storage (allocates one extra field).
    #[must_use]
    pub fn with_pressure(mut self) -> Self {
        self.pressure = Some(vec![T::zero(); self.nx * self.ny]);
        self
    }

    /// Enable nuclei fraction tracking (allocates one extra field).
    #[must_use]
    pub fn with_nuclei(mut self) -> Self {
        self.nuclei_fraction = Some(vec![T::zero(); self.nx * self.ny]);
        self
    }

    /// Update all macroscopic fields from the flat distribution buffer.
    ///
    /// # Invariant (Theorem — Moment Consistency)
    /// After this call, `density[j, i]` equals $\sum_q f_q(i,j)$ and
    /// `velocity[j, i]` equals $(\sum_q e_{q,x} f_q) / \rho$, both exact.
    pub fn update_from_distributions(&mut self, f: &[T], g: Option<&[T]>) {
        let nx = self.nx;
        let ny = self.ny;

        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;
                let rho = compute_density_flat(f, j, i, nx);
                self.density[cell] = rho;
                let [ux, uy] = compute_velocity_flat::<T>(f, j, i, nx, rho);
                self.velocity[cell * 2] = ux;
                self.velocity[cell * 2 + 1] = uy;

                if let Some(ref mut pressure) = self.pressure {
                    pressure[cell] = compute_pressure(rho);
                }

                if let (Some(g_slice), Some(ref mut nuclei)) = (g, &mut self.nuclei_fraction) {
                    nuclei[cell] = compute_density_flat(g_slice, j, i, nx);
                }
            }
        }
    }

    /// Get density at node (i, j).
    #[inline]
    pub fn density_at(&self, i: usize, j: usize) -> T {
        self.density[j * self.nx + i]
    }

    /// Get velocity at node (i, j) as [u_x, u_y].
    #[inline]
    pub fn velocity_at(&self, i: usize, j: usize) -> [T; 2] {
        let base = (j * self.nx + i) * 2;
        [self.velocity[base], self.velocity[base + 1]]
    }
}

/// Compute density at node (i, j) from the flat buffer.
///
/// ρ(i,j) = ∑_{q=0}^{8} f_q(i,j)
#[inline]
pub fn compute_density_flat<T: RealField + Copy>(f: &[T], j: usize, i: usize, nx: usize) -> T {
    let mut rho = T::zero();
    for q in 0..9 {
        rho += f[f_idx(j, i, q, nx)];
    }
    rho
}

/// Compute velocity at node (i, j) from the flat buffer.
///
/// **u**(i,j) = (∑_q **e**_q f_q(i,j)) / ρ(i,j)
#[inline]
pub fn compute_velocity_flat<T: RealField + Copy + FromPrimitive>(
    f: &[T],
    j: usize,
    i: usize,
    nx: usize,
    density: T,
) -> [T; 2] {
    let mut ux = T::zero();
    let mut uy = T::zero();
    for q in 0..9 {
        let (ex, ey) = D2Q9::VELOCITIES[q];
        let fq = f[f_idx(j, i, q, nx)];
        ux += T::from_i32(ex).expect("lattice velocity is a small integer") * fq;
        uy += T::from_i32(ey).expect("lattice velocity is a small integer") * fq;
    }
    [ux / density, uy / density]
}

/// Compute density from a node-local 9-element slice.
///
/// Used for backward-compatible helper contexts where the 9 values are
/// already extracted into a `[T; 9]`.
pub fn compute_density<T: RealField + Copy>(f_node: &[T; 9]) -> T {
    f_node.iter().fold(T::zero(), |acc, &v| acc + v)
}

/// Compute velocity from a node-local slice (for helper contexts).
pub fn compute_velocity<T: RealField + Copy + FromPrimitive>(
    f_node: &[T; 9],
    density: T,
) -> [T; 2] {
    let mut ux = T::zero();
    let mut uy = T::zero();
    for q in 0..9 {
        let (ex, ey) = D2Q9::VELOCITIES[q];
        let ex_t = T::from_i32(ex).expect("lattice velocity is a small integer");
        let ey_t = T::from_i32(ey).expect("lattice velocity is a small integer");
        ux += ex_t * f_node[q];
        uy += ey_t * f_node[q];
    }
    [ux / density, uy / density]
}

/// LBM equation of state: p = c_s² ρ with c_s² = 1/3 (Theorem — LBM Pressure Relation).
pub fn compute_pressure<T: RealField + Copy + FromPrimitive>(density: T) -> T {
    let cs2 = T::from_f64(1.0 / 3.0)
        .expect("T must represent f64 values; cs² = 1/3 is exact in IEEE 754");
    cs2 * density
}

/// Compute stress tensor Π_{αβ} = ∑_q e_{q,α} e_{q,β} (f_q - f_q^eq).
pub fn compute_stress_tensor<T: RealField + Copy + FromPrimitive>(
    f_node: &[T; 9],
    f_eq: &[T; 9],
) -> [[T; 2]; 2] {
    let mut stress = [[T::zero(); 2]; 2];
    for q in 0..9 {
        let (ex, ey) = D2Q9::VELOCITIES[q];
        let ex_t = T::from_i32(ex).expect("lattice velocity is a small integer");
        let ey_t = T::from_i32(ey).expect("lattice velocity is a small integer");
        let f_neq = f_node[q] - f_eq[q];
        stress[0][0] += ex_t * ex_t * f_neq;
        stress[0][1] += ex_t * ey_t * f_neq;
        stress[1][0] += ey_t * ex_t * f_neq;
        stress[1][1] += ey_t * ey_t * f_neq;
    }
    stress
}

/// Compute kinetic energy density: E_k = ½ ρ |**u**|².
pub fn compute_kinetic_energy<T: RealField + Copy + FromPrimitive>(
    density: T,
    velocity: &[T; 2],
) -> T {
    let half = T::from_f64(0.5).expect("T must represent f64 values; ½ is exact in IEEE 754");
    half * density * (velocity[0] * velocity[0] + velocity[1] * velocity[1])
}

/// Compute vorticity ω = ∂v/∂x − ∂u/∂y at node (i, j) using central differences.
pub fn compute_vorticity<T: RealField + Copy + FromPrimitive>(
    velocity: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    i: usize,
    j: usize,
) -> T {
    if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
        return T::zero();
    }
    let two = T::from_f64(2.0).expect("T must represent f64 values");
    let cell = |jj: usize, ii: usize, d: usize| velocity[(jj * nx + ii) * 2 + d];
    let dv_dx = (cell(j, i + 1, 1) - cell(j, i - 1, 1)) / (two * dx);
    let du_dy = (cell(j + 1, i, 0) - cell(j - 1, i, 0)) / (two * dy);
    dv_dx - du_dy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_computation() {
        let f_node = [0.1_f64; 9];
        let density = compute_density(&f_node);
        assert!((density - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_computation() {
        let mut f_node = [0.1_f64; 9];
        f_node[1] = 0.15; // East (ex=+1)
        f_node[3] = 0.05; // West (ex=-1)
        let density = compute_density(&f_node);
        let velocity = compute_velocity(&f_node, density);
        assert!(velocity[0] > 0.0);
        assert!(velocity[1].abs() < 1e-10);
    }

    #[test]
    fn test_pressure_computation() {
        let density = 1.5_f64;
        let pressure = compute_pressure(density);
        let expected = (1.0 / 3.0) * density;
        assert!((pressure - expected).abs() < 1e-10);
    }

    #[test]
    fn test_flat_density_matches_node_density() {
        let nx = 4_usize;
        let ny = 4_usize;
        let mut f = vec![0.0_f64; nx * ny * 9];
        // Set node (i=2, j=1) to non-trivial values
        for q in 0..9 {
            f[f_idx(1, 2, q, nx)] = (q as f64 + 1.0) * 0.05;
        }
        let rho_flat = compute_density_flat(&f, 1, 2, nx);
        let rho_node_sum: f64 = (0..9).map(|q| (q as f64 + 1.0) * 0.05).sum();
        assert!((rho_flat - rho_node_sum).abs() < 1e-12);
    }
}
