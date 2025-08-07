//! Conservation checking for CFD simulations.
//!
//! This module provides tools to verify that CFD simulations satisfy fundamental
//! conservation laws such as mass, momentum, and energy conservation.

use cfd_core::Result;
use nalgebra::{RealField, DVector};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Trait for conservation checking
pub trait ConservationChecker<T: RealField> {
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
pub struct ConservationReport<T: RealField> {
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
}

/// Tolerance settings for conservation checks
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationTolerance<T: RealField> {
    /// Absolute tolerance
    pub absolute: T,
    /// Relative tolerance
    pub relative: T,
}

impl<T: RealField + FromPrimitive> Default for ConservationTolerance<T> {
    fn default() -> Self {
        Self {
            absolute: T::from_f64(1e-12).unwrap(),
            relative: T::from_f64(1e-10).unwrap(),
        }
    }
}

/// History of conservation errors over time
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationHistory<T: RealField> {
    /// Time points
    pub times: Vec<T>,
    /// Conservation errors at each time point
    pub errors: Vec<T>,
    /// Whether conservation was satisfied at each time point
    pub satisfied: Vec<bool>,
}

impl<T: RealField> ConservationHistory<T> {
    /// Create a new conservation history
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            errors: Vec::new(),
            satisfied: Vec::new(),
        }
    }

    /// Add a new time point to the history
    pub fn add_point(&mut self, time: T, error: T, satisfied: bool) {
        self.times.push(time);
        self.errors.push(error);
        self.satisfied.push(satisfied);
    }

    /// Get the maximum error in the history
    pub fn max_error(&self) -> Option<T> {
        self.errors.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).cloned()
    }

    /// Get the fraction of time points where conservation was satisfied
    pub fn satisfaction_rate(&self) -> f64 {
        if self.satisfied.is_empty() {
            0.0
        } else {
            let satisfied_count = self.satisfied.iter().filter(|&&x| x).count();
            satisfied_count as f64 / self.satisfied.len() as f64
        }
    }
}

/// Mass conservation checker
#[derive(Debug)]
pub struct MassConservation<T: RealField> {
    /// Tolerance for mass conservation
    tolerance: ConservationTolerance<T>,
}

impl<T: RealField + FromPrimitive> MassConservation<T> {
    /// Create a new mass conservation checker
    pub fn new(tolerance: ConservationTolerance<T>) -> Self {
        Self { tolerance }
    }

    /// Create with default tolerance
    pub fn default() -> Self {
        Self::new(ConservationTolerance::default())
    }
}

impl<T: RealField + FromPrimitive> ConservationChecker<T> for MassConservation<T> {
    type FlowField = (DVector<T>, DVector<T>); // (density, velocity_divergence)

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (density, velocity_div) = field;

        // Check continuity equation: ∂ρ/∂t + ∇·(ρv) = 0
        // For steady flow: ∇·(ρv) = 0
        // Simplified check: sum of mass flux should be zero
        let mass_flux_sum = density.iter()
            .zip(velocity_div.iter())
            .map(|(rho, div_v)| rho.clone() * div_v.clone())
            .fold(T::zero(), |acc, x| acc + x);

        let error = mass_flux_sum.clone().abs();
        let is_conserved = error <= self.tolerance.absolute;

        let mut details = HashMap::new();
        details.insert("mass_flux_sum".to_string(), mass_flux_sum);
        details.insert("absolute_error".to_string(), error.clone());

        Ok(ConservationReport {
            check_name: self.name().to_string(),
            is_conserved,
            error,
            tolerance: self.tolerance.absolute.clone(),
            details,
        })
    }

    fn name(&self) -> &str {
        "Mass Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance.absolute.clone()
    }
}

/// Energy conservation checker
#[derive(Debug)]
pub struct EnergyConservation<T: RealField> {
    /// Tolerance for energy conservation
    tolerance: ConservationTolerance<T>,
}

impl<T: RealField + FromPrimitive> EnergyConservation<T> {
    /// Create a new energy conservation checker
    pub fn new(tolerance: ConservationTolerance<T>) -> Self {
        Self { tolerance }
    }

    /// Create with default tolerance
    pub fn default() -> Self {
        Self::new(ConservationTolerance::default())
    }
}

impl<T: RealField + FromPrimitive + std::iter::Sum> ConservationChecker<T> for EnergyConservation<T> {
    type FlowField = (DVector<T>, DVector<T>, DVector<T>); // (kinetic_energy, potential_energy, dissipation)

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (kinetic, potential, dissipation) = field;

        // Simplified energy balance check
        let total_energy: T = kinetic.iter().zip(potential.iter()).map(|(k, p)| k.clone() + p.clone()).sum();
        let total_dissipation: T = dissipation.iter().cloned().sum();

        // For steady flow, energy input should equal dissipation
        let energy_imbalance = (total_energy.clone() - total_dissipation.clone()).abs();
        let relative_error = if total_energy > T::zero() {
            energy_imbalance.clone() / total_energy.clone()
        } else {
            energy_imbalance.clone()
        };

        let is_conserved = energy_imbalance <= self.tolerance.absolute
            && relative_error <= self.tolerance.relative;

        let mut details = HashMap::new();
        details.insert("total_energy".to_string(), total_energy);
        details.insert("total_dissipation".to_string(), total_dissipation);
        details.insert("energy_imbalance".to_string(), energy_imbalance.clone());
        details.insert("relative_error".to_string(), relative_error);

        Ok(ConservationReport {
            check_name: self.name().to_string(),
            is_conserved,
            error: energy_imbalance,
            tolerance: self.tolerance.absolute.clone(),
            details,
        })
    }

    fn name(&self) -> &str {
        "Energy Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance.absolute.clone()
    }
}
