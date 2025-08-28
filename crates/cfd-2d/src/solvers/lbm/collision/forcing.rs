//! Forcing schemes for external forces in LBM

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Forcing scheme trait
pub trait ForcingScheme<T: RealField + Copy> {
    /// Apply forcing to distribution functions
    fn apply_force(&self, f: &mut [T; 9], force: [T; 2], velocity: [T; 2], density: T, dt: T);
}

/// Guo forcing scheme
/// Based on Guo et al. (2002)
pub struct GuoForcing<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> GuoForcing<T> {
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ForcingScheme<T> for GuoForcing<T> {
    fn apply_force(&self, f: &mut [T; 9], force: [T; 2], velocity: [T; 2], density: T, dt: T) {
        use crate::solvers::lbm::lattice::D2Q9;

        let cs2 =
            T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));

        for q in 0..9 {
            let c = D2Q9::velocity_vector::<T>(q);
            let w = D2Q9::weight::<T>(q);

            // Guo forcing term
            let cu = c[0] * velocity[0] + c[1] * velocity[1];
            let cf = c[0] * force[0] + c[1] * force[1];

            let forcing_term = w * density * dt * (cf / cs2 + (cu * cf) / (cs2 * cs2));

            f[q] = f[q] + forcing_term;
        }
    }
}

/// Shan-Chen forcing scheme
/// For multiphase flows
pub struct ShanChenForcing<T: RealField + Copy> {
    /// Interaction strength
    g: T,
}

impl<T: RealField + Copy> ShanChenForcing<T> {
    pub fn new(interaction_strength: T) -> Self {
        Self {
            g: interaction_strength,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ForcingScheme<T> for ShanChenForcing<T> {
    fn apply_force(&self, f: &mut [T; 9], force: [T; 2], velocity: [T; 2], density: T, dt: T) {
        use crate::solvers::lbm::lattice::D2Q9;

        // Shan-Chen forcing implementation
        for q in 0..9 {
            let c = D2Q9::velocity_vector::<T>(q);
            let w = D2Q9::weight::<T>(q);

            // Interaction force
            let cf = c[0] * force[0] + c[1] * force[1];
            let forcing_term = self.g * w * cf * dt;

            f[q] = f[q] + forcing_term;
        }
    }
}
