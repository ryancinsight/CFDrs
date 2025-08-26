//! Conservation checking for CFD simulations.
//!
//! This module provides tools to verify that CFD simulations satisfy fundamental
//! conservation laws such as mass, momentum, and energy conservation.

use cfd_core::error::Result;
use cfd_core::numeric;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
/// Trait for conservation checking
pub trait ConservationChecker<T: RealField + Copy> {
    /// Type representing the flow field
    type FlowField;
    /// Check conservation for the given flow field
    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>>;
    /// Get the name of this conservation check
    fn name(&self) -> &str;
    /// Get the tolerance for this conservation check
    fn tolerance(&self) -> T;
}
/// Report on conservation properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationReport<T: RealField + Copy> {
    /// Name of the conservation check
    pub check_name: String,
    /// Whether conservation is satisfied within tolerance
    pub is_conserved: bool,
    /// Conservation error magnitude
    pub error: T,
    /// Tolerance used for the check
    pub tolerance: T,
    /// Additional details about the conservation check
    pub details: HashMap<String, T>,
/// Tolerance settings for conservation checks
    }

}

pub struct ConservationTolerance<T: RealField + Copy> {
    /// Absolute tolerance
    pub absolute: T,
    /// Relative tolerance
    pub relative: T,
impl<T: RealField + Copy + FromPrimitive + Copy> Default for ConservationTolerance<T> {
    fn default() -> Self {
        Self {
            absolute: cfd_core::numeric::from_f64(1e-12)?,
            relative: cfd_core::numeric::from_f64(1e-10)?,
        }
    }
/// History of conservation errors over time
pub struct ConservationHistory<T: RealField + Copy> {
    /// Time points
    pub times: Vec<T>,
    /// Conservation errors at each time point
    pub errors: Vec<T>,
    /// Whether conservation was satisfied at each time point
    pub satisfied: Vec<bool>,
impl<T: RealField + Copy> Default for ConservationHistory<T> {
        Self::new()
impl<T: RealField + Copy> ConservationHistory<T> {
    /// Create a new conservation history
    #[must_use]
    pub fn new() -> Self {
            times: Vec::new(),
            errors: Vec::new(),
            satisfied: Vec::new(),
    /// Add a new time point to the history
    }

    pub fn add_point(&mut self, time: T, error: T, satisfied: bool) {
        self.times.push(time);
        self.errors.push(error);
        self.satisfied.push(satisfied);
    /// Get the maximum error in the history
    }

    pub fn max_error(&self) -> Option<T> {
        self.errors
            .iter()
            .max_by(|a, b| {
                a.partial_cmp(b)
                    .expect("CRITICAL: Add proper error handling")
            })
            .copied()
    /// Get the fraction of time points where conservation was satisfied
    pub fn satisfaction_rate(&self) -> f64 {
        if self.satisfied.is_empty() {
            0.0
        } else {
            let satisfied_count = self.satisfied.iter().filter(|&&x| x).count();
            satisfied_count as f64 / self.satisfied.len() as f64
/// Mass conservation checker
#[derive(Debug)]
pub struct MassConservation<T: RealField + Copy> {
    /// Tolerance for mass conservation
    tolerance: ConservationTolerance<T>,
impl<T: RealField + Copy + FromPrimitive + Copy> MassConservation<T> {
    /// Create a new mass conservation checker
    pub fn new(tolerance: ConservationTolerance<T>) -> Self {
        Self { tolerance }
    /// Create with default tolerance
    pub fn default() -> Self {
        Self::new(ConservationTolerance::default())
impl<T: RealField + Copy + FromPrimitive + Copy> ConservationChecker<T> for MassConservation<T> {
    type FlowField = (DVector<T>, DVector<T>); // (density, velocity_divergence)
    }

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (density, velocity_div) = field;
        // Check continuity equation: ∂ρ/∂t + ∇·(ρv) = 0
        // For steady incompressible flow: ∇·v = 0
        // For compressible flow: ∇·(ρv) = 0
        let mass_flux_sum = density
            .zip(velocity_div.iter())
            .map(|(rho, div_v)| *rho * *div_v)
            .fold(T::zero(), |acc, x| acc + x);
        let error = mass_flux_sum.abs();
        let is_conserved = error <= self.tolerance.absolute;
        let mut details = HashMap::new();
        details.insert("mass_flux_sum".to_string(), mass_flux_sum);
        details.insert("absolute_error".to_string(), error);
        Ok(ConservationReport {
            check_name: self.name().to_string(),
            is_conserved,
            error,
            tolerance: self.tolerance.absolute,
            details,
        })
    fn name(&self) -> &str {
        "Mass Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance.absolute
/// Energy conservation checker
    }

pub struct EnergyConservation<T: RealField + Copy> {
    /// Tolerance for energy conservation
impl<T: RealField + Copy + FromPrimitive + Copy> EnergyConservation<T> {
    /// Create a new energy conservation checker
impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> ConservationChecker<T>
    for EnergyConservation<T>
{
    type FlowField = (DVector<T>, DVector<T>, DVector<T>); // (kinetic_energy, potential_energy, dissipation)
        let (kinetic, potential, dissipation) = field;
        // Energy balance check: total energy = kinetic + potential
        let total_energy: T = kinetic
            .zip(potential.iter())
            .map(|(k, p)| *k + *p)
            .sum();
        let total_dissipation: T = dissipation.iter().copied().sum();
        // For steady flow, energy input should equal dissipation
        let energy_imbalance = (total_energy - total_dissipation).abs();
        let relative_error = if total_energy > T::zero() {
            energy_imbalance / total_energy
            energy_imbalance
        };
        let is_conserved = energy_imbalance <= self.tolerance.absolute
            && relative_error <= self.tolerance.relative;
        details.insert("total_energy".to_string(), total_energy);
        details.insert("total_dissipation".to_string(), total_dissipation);
        details.insert("energy_imbalance".to_string(), energy_imbalance);
        details.insert("relative_error".to_string(), relative_error);
            error: energy_imbalance,
        "Energy Conservation"
/// Global conservation integrals for comprehensive validation
#[derive(Debug, Clone)]
pub struct GlobalConservationIntegrals<T: RealField + Copy> {
    /// Total mass in the domain
    pub total_mass: T,
    /// Total momentum in each direction
    pub total_momentum: [T; 3],
    /// Total kinetic energy
    pub total_kinetic_energy: T,
    /// Domain volume/area
    pub domain_volume: T,
    /// Mass flux through boundaries
    pub boundary_mass_flux: T,
    /// Momentum flux through boundaries
    pub boundary_momentum_flux: [T; 3],
    /// Energy flux through boundaries
    pub boundary_energy_flux: T,
impl<T: RealField + Copy + FromPrimitive + Copy> Default for GlobalConservationIntegrals<T> {
impl<T: RealField + Copy + FromPrimitive + Copy> GlobalConservationIntegrals<T> {
    /// Create new global conservation integrals
            total_mass: T::zero(),
            total_momentum: [T::zero(), T::zero(), T::zero()],
            total_kinetic_energy: T::zero(),
            domain_volume: T::zero(),
            boundary_mass_flux: T::zero(),
            boundary_momentum_flux: [T::zero(), T::zero(), T::zero()],
            boundary_energy_flux: T::zero(),
    /// Compute conservation integrals from velocity and density fields using iterator optimization
    pub fn compute_from_fields(
        density: &[T],
        velocity_x: &[T],
        velocity_y: &[T],
        velocity_z: &[T],
        cell_volumes: &[T],
    ) -> Result<Self> {
        if density.len() != velocity_x.len()
            || density.len() != velocity_y.len()
            || density.len() != velocity_z.len()
            || density.len() != cell_volumes.len()
        {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "All field arrays must have the same length".to_string(),
            ));
        let mut integrals = Self::new();
        // Use iterator combinators for zero-copy computation
        let (total_mass, total_momentum, total_kinetic_energy, domain_volume) = density
            .zip(velocity_x.iter())
            .zip(velocity_y.iter())
            .zip(velocity_z.iter())
            .zip(cell_volumes.iter())
            .fold(
                (
                    T::zero(),
                    [T::zero(), T::zero(), T::zero()],
                ),
                |(mass, [mx, my, mz], ke, vol), ((((rho, u), v), w), dv)| {
                    let dm = *rho * *dv;
                    let u_sq = *u * *u;
                    let v_sq = *v * *v;
                    let w_sq = *w * *w;
                    let half = cfd_core::numeric::from_f64(0.5)?;
                    (
                        mass + dm,
                        [mx + dm * *u, my + dm * *v, mz + dm * *w],
                        ke + dm * (u_sq + v_sq + w_sq) * half,
                        vol + *dv,
                    )
                },
            );
        integrals.total_mass = total_mass;
        integrals.total_momentum = total_momentum;
        integrals.total_kinetic_energy = total_kinetic_energy;
        integrals.domain_volume = domain_volume;
        Ok(integrals)
    /// Check mass conservation: dm/dt + ∇·(ρu) = 0
    pub fn mass_conservation_error(&self, previous: &Self, dt: T) -> T {
        let mass_change = (self.total_mass - previous.total_mass) / dt;
        (mass_change + self.boundary_mass_flux).abs()
    /// Check momentum conservation: d(ρu)/dt + ∇·(ρuu) = -∇p + ∇·τ + F
    }

    pub fn momentum_conservation_error(&self, previous: &Self, dt: T) -> [T; 3] {
        [
            ((self.total_momentum[0] - previous.total_momentum[0]) / dt
                + self.boundary_momentum_flux[0])
                .abs(),
            ((self.total_momentum[1] - previous.total_momentum[1]) / dt
                + self.boundary_momentum_flux[1])
            ((self.total_momentum[2] - previous.total_momentum[2]) / dt
                + self.boundary_momentum_flux[2])
        ]
    /// Check energy conservation: dE/dt + ∇·(u(E+p)) = ∇·(k∇T) + Φ
    }

    pub fn energy_conservation_error(&self, previous: &Self, dt: T) -> T {
        let energy_change = (self.total_kinetic_energy - previous.total_kinetic_energy) / dt;
        (energy_change + self.boundary_energy_flux).abs()

    }


}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
