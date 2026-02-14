//! Vascular bifurcation network models
//!
//! This module provides 1D network models for arterial and venous bifurcations,
//! supporting hierarchical vascular trees with pressure-flow coupling.
//!
//! # Mathematical Foundation
//!
//! ## Junction Conditions
//! At a bifurcation junction:
//! 1. Mass conservation: Q₀ = Q₁ + Q₂
//! 2. Pressure continuity: p₀ = p₁ = p₂ (no junction loss)
//!    OR p₁ = p₀ - K·ρ·V²/2 (with junction loss)
//!
//! ## Network Resistance
//! Series: R_total = R₁ + R₂
//! Parallel: 1/R_total = 1/R₁ + 1/R₂
//!
//! # References
//! - Olufsen, M.S. (1999) "Structured tree outflow condition"
//! - Sherwin, S.J. et al. (2003) "One-dimensional modelling of a vascular network"

use super::murrays_law::MurraysLaw;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// Junction Types
// ============================================================================

/// Type of vascular junction
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionType {
    /// Simple bifurcation (one parent, two daughters)
    Bifurcation,
    /// Trifurcation (one parent, three daughters)
    Trifurcation,
    /// Confluence (two parents, one daughter) - venous junction
    Confluence,
    /// Anastomosis (connecting collateral)
    Anastomosis,
}

/// Junction loss model
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionLossModel {
    /// No pressure loss at junction (ideal)
    None,
    /// Simple K-factor loss: ΔP = K · ρV²/2
    KFactor,
    /// Energy-preserving junction
    EnergyPreserving,
}

// ============================================================================
// Vessel Segment
// ============================================================================

/// A single vessel segment in the vascular network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VesselSegment<T: RealField + Copy> {
    /// Segment identifier
    pub id: usize,
    /// Vessel radius [m]
    pub radius: T,
    /// Vessel length [m]
    pub length: T,
    /// Wall thickness [m]
    pub wall_thickness: T,
    /// Young's modulus of wall [Pa]
    pub youngs_modulus: T,
    /// Inlet node ID
    pub inlet_node: usize,
    /// Outlet node ID
    pub outlet_node: usize,
}

impl<T: RealField + FromPrimitive + Copy> VesselSegment<T> {
    /// Create new vessel segment
    pub fn new(
        id: usize,
        radius: T,
        length: T,
        inlet_node: usize,
        outlet_node: usize,
    ) -> Self {
        Self {
            id,
            radius,
            length,
            wall_thickness: radius * T::from_f64(0.1).unwrap(), // 10% of radius
            youngs_modulus: T::from_f64(0.4e6).unwrap(),        // ~0.4 MPa for arteries
            inlet_node,
            outlet_node,
        }
    }

    /// Calculate Poiseuille resistance R = 8μL/(πR⁴)
    pub fn resistance(&self, viscosity: T) -> T {
        let pi = T::from_f64(PI).unwrap();
        let eight = T::from_f64(8.0).unwrap();
        eight * viscosity * self.length / (pi * self.radius.powi(4))
    }

    /// Calculate inertance L = ρL/(πR²)
    pub fn inertance(&self, density: T) -> T {
        let pi = T::from_f64(PI).unwrap();
        density * self.length / (pi * self.radius * self.radius)
    }

    /// Calculate compliance C = 3πR³L/(2Eh)
    /// Based on thin-walled tube approximation
    pub fn compliance(&self) -> T {
        let pi = T::from_f64(PI).unwrap();
        let three = T::from_f64(3.0).unwrap();
        let two = T::from_f64(2.0).unwrap();
        three * pi * self.radius.powi(3) * self.length
            / (two * self.youngs_modulus * self.wall_thickness)
    }

    /// Calculate wave speed c = √(Eh/(2ρR))
    pub fn wave_speed(&self, density: T) -> T {
        let two = T::from_f64(2.0).unwrap();
        (self.youngs_modulus * self.wall_thickness / (two * density * self.radius)).sqrt()
    }

    /// Calculate cross-sectional area
    pub fn area(&self) -> T {
        let pi = T::from_f64(PI).unwrap();
        pi * self.radius * self.radius
    }

    /// Calculate diameter
    pub fn diameter(&self) -> T {
        let two = T::from_f64(2.0).unwrap();
        two * self.radius
    }
}

// ============================================================================
// Bifurcation Node
// ============================================================================

/// A bifurcation junction connecting vessel segments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bifurcation<T: RealField + Copy> {
    /// Node identifier
    pub id: usize,
    /// Junction type
    pub junction_type: JunctionType,
    /// Loss model at junction
    pub loss_model: JunctionLossModel,
    /// K-factor for pressure loss (if applicable)
    pub k_factor: T,
    /// Incoming vessel IDs
    pub parent_vessels: Vec<usize>,
    /// Outgoing vessel IDs
    pub daughter_vessels: Vec<usize>,
    /// Current pressure at junction [Pa]
    pub pressure: T,
    /// Current flow rates [m³/s] - indexed by vessel ID
    pub flow_rates: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> Bifurcation<T> {
    /// Create new bifurcation junction
    pub fn new(id: usize, parent_id: usize, daughter_ids: Vec<usize>) -> Self {
        let n_flows = 1 + daughter_ids.len();
        Self {
            id,
            junction_type: if daughter_ids.len() == 2 {
                JunctionType::Bifurcation
            } else {
                JunctionType::Trifurcation
            },
            loss_model: JunctionLossModel::None,
            k_factor: T::from_f64(0.1).unwrap(),
            parent_vessels: vec![parent_id],
            daughter_vessels: daughter_ids,
            pressure: T::from_f64(13_332.0).unwrap(), // ~100 mmHg
            flow_rates: vec![T::zero(); n_flows],
        }
    }

    /// Create confluence (venous) junction
    pub fn confluence(id: usize, parent_ids: Vec<usize>, daughter_id: usize) -> Self {
        let n_flows = parent_ids.len() + 1;
        Self {
            id,
            junction_type: JunctionType::Confluence,
            loss_model: JunctionLossModel::None,
            k_factor: T::from_f64(0.05).unwrap(),
            parent_vessels: parent_ids,
            daughter_vessels: vec![daughter_id],
            pressure: T::from_f64(1_333.0).unwrap(), // ~10 mmHg for veins
            flow_rates: vec![T::zero(); n_flows],
        }
    }

    /// Calculate mass conservation error
    pub fn mass_conservation_error(&self) -> T {
        if self.flow_rates.is_empty() {
            return T::zero();
        }

        let parent_flow: T = self.flow_rates.iter().take(self.parent_vessels.len()).fold(T::zero(), |acc, &f| acc + f);
        let daughter_flow: T = self.flow_rates.iter().skip(self.parent_vessels.len()).fold(T::zero(), |acc, &f| acc + f);

        if parent_flow.abs() < T::from_f64(1e-20).unwrap() {
            return T::zero();
        }

        (parent_flow - daughter_flow).abs() / parent_flow.abs()
    }

    /// Apply junction loss model to get pressure at daughters
    pub fn daughter_pressure(&self, parent_pressure: T, parent_velocity: T, density: T) -> T {
        match self.loss_model {
            JunctionLossModel::None => parent_pressure,
            JunctionLossModel::KFactor => {
                let half = T::from_f64(0.5).unwrap();
                parent_pressure - self.k_factor * half * density * parent_velocity * parent_velocity
            }
            JunctionLossModel::EnergyPreserving => {
                // Energy preservation: p + ρv²/2 = const
                // More complex model would account for velocity changes
                parent_pressure
            }
        }
    }
}

// ============================================================================
// Bifurcation Network
// ============================================================================

/// A network of vessel segments connected by bifurcations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationNetwork<T: RealField + Copy> {
    /// All vessel segments
    pub vessels: Vec<VesselSegment<T>>,
    /// All junction nodes
    pub junctions: Vec<Bifurcation<T>>,
    /// Inlet boundary pressure [Pa]
    pub inlet_pressure: T,
    /// Outlet boundary resistance [Pa·s/m³]
    pub outlet_resistance: T,
    /// Blood density [kg/m³]
    pub density: T,
    /// Blood viscosity [Pa·s]
    pub viscosity: T,
}

impl<T: RealField + FromPrimitive + Copy> BifurcationNetwork<T> {
    /// Create empty network with default blood properties
    pub fn new() -> Self {
        Self {
            vessels: Vec::new(),
            junctions: Vec::new(),
            inlet_pressure: T::from_f64(13_332.0).unwrap(), // 100 mmHg
            outlet_resistance: T::from_f64(1e9).unwrap(),
            density: T::from_f64(1060.0).unwrap(),
            viscosity: T::from_f64(0.0035).unwrap(),
        }
    }

    /// Add a vessel segment
    pub fn add_vessel(&mut self, vessel: VesselSegment<T>) -> usize {
        let id = self.vessels.len();
        self.vessels.push(vessel);
        id
    }

    /// Add a bifurcation junction
    pub fn add_junction(&mut self, junction: Bifurcation<T>) -> usize {
        let id = self.junctions.len();
        self.junctions.push(junction);
        id
    }

    /// Create a Murray's Law compliant symmetric bifurcation tree
    ///
    /// # Arguments
    /// * `root_radius` - Root vessel radius [m]
    /// * `root_length` - Root vessel length [m]
    /// * `generations` - Number of bifurcation generations (1 = just root)
    /// * `length_ratio` - Length of daughter / length of parent
    pub fn create_symmetric_tree(
        root_radius: T,
        root_length: T,
        generations: usize,
        length_ratio: T,
    ) -> Self {
        let mut network = Self::new();
        let murray = MurraysLaw::<T>::new();

        // Recursive helper to build tree
        fn add_generation<T: RealField + FromPrimitive + Copy>(
            network: &mut BifurcationNetwork<T>,
            murray: &MurraysLaw<T>,
            _parent_id: usize,
            parent_radius: T,
            parent_length: T,
            length_ratio: T,
            gen: usize,
            max_gen: usize,
            inlet_node: usize,
        ) -> usize {
            let outlet_node = network.vessels.len() + 1;

            // Add parent vessel
            let vessel = VesselSegment::new(
                network.vessels.len(),
                parent_radius,
                parent_length,
                inlet_node,
                outlet_node,
            );
            network.vessels.push(vessel);

            if gen >= max_gen {
                return outlet_node;
            }

            // Create daughters
            let daughter_radius = murray.symmetric_daughter_diameter(parent_radius * T::from_f64(2.0).unwrap()) / T::from_f64(2.0).unwrap();
            let daughter_length = parent_length * length_ratio;

            let parent_vessel_id = network.vessels.len() - 1;

            // Add junction
            let daughter1_id = network.vessels.len();

            // Recursively add daughters
            let _terminal1 = add_generation(
                network,
                murray,
                daughter1_id,
                daughter_radius,
                daughter_length,
                length_ratio,
                gen + 1,
                max_gen,
                outlet_node,
            );

            let daughter2_actual_id = network.vessels.len();
            let _terminal2 = add_generation(
                network,
                murray,
                daughter2_actual_id,
                daughter_radius,
                daughter_length,
                length_ratio,
                gen + 1,
                max_gen,
                outlet_node,
            );

            // Create bifurcation junction
            let junction = Bifurcation::new(
                outlet_node,
                parent_vessel_id,
                vec![daughter1_id, daughter2_actual_id],
            );
            network.junctions.push(junction);

            outlet_node
        }

        add_generation(
            &mut network,
            &murray,
            0,
            root_radius,
            root_length,
            length_ratio,
            1,
            generations,
            0,
        );

        network
    }

    /// Calculate total network resistance
    pub fn total_resistance(&self) -> T {
        // Simple sum of all vessel resistances (series approximation)
        self.vessels
            .iter()
            .map(|v| v.resistance(self.viscosity))
            .fold(T::zero(), |acc, r| acc + r)
    }

    /// Get number of terminal vessels
    pub fn terminal_count(&self) -> usize {
        // Vessels that don't appear as parents in any junction
        let parent_vessels: std::collections::HashSet<_> = self
            .junctions
            .iter()
            .flat_map(|j| j.parent_vessels.iter())
            .collect();

        self.vessels
            .iter()
            .filter(|v| !parent_vessels.contains(&v.id))
            .count()
    }

    /// Validate Murray's Law compliance for all junctions
    pub fn validate_murrays_law(&self, tolerance: T) -> bool {
        let murray = MurraysLaw::<T>::new();

        for junction in &self.junctions {
            if junction.junction_type != JunctionType::Bifurcation {
                continue;
            }

            if junction.parent_vessels.len() != 1 || junction.daughter_vessels.len() != 2 {
                continue;
            }

            let parent = &self.vessels[junction.parent_vessels[0]];
            let d1 = &self.vessels[junction.daughter_vessels[0]];
            let d2 = &self.vessels[junction.daughter_vessels[1]];

            if !murray.is_valid(parent.diameter(), d1.diameter(), d2.diameter(), tolerance) {
                return false;
            }
        }

        true
    }
}

impl<T: RealField + FromPrimitive + Copy> Default for BifurcationNetwork<T> {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vessel_segment_resistance() {
        let vessel = VesselSegment::<f64>::new(0, 0.003, 0.1, 0, 1);
        let r = vessel.resistance(0.0035);

        // R = 8μL/(πR⁴) should be finite and positive
        assert!(r > 0.0 && r.is_finite());
    }

    #[test]
    fn test_vessel_wave_speed() {
        let vessel = VesselSegment::<f64>::new(0, 0.003, 0.1, 0, 1);
        let c = vessel.wave_speed(1060.0);

        // Wave speed should be ~5-10 m/s for arteries
        assert!(c > 1.0 && c < 20.0, "Wave speed {} should be 5-10 m/s", c);
    }

    #[test]
    fn test_bifurcation_creation() {
        let bif = Bifurcation::<f64>::new(0, 0, vec![1, 2]);

        assert_eq!(bif.junction_type, JunctionType::Bifurcation);
        assert_eq!(bif.parent_vessels.len(), 1);
        assert_eq!(bif.daughter_vessels.len(), 2);
    }

    #[test]
    fn test_confluence_creation() {
        let conf = Bifurcation::<f64>::confluence(0, vec![0, 1], 2);

        assert_eq!(conf.junction_type, JunctionType::Confluence);
        assert_eq!(conf.parent_vessels.len(), 2);
        assert_eq!(conf.daughter_vessels.len(), 1);
    }

    #[test]
    fn test_network_creation() {
        let network = BifurcationNetwork::<f64>::new();

        assert!(network.vessels.is_empty());
        assert!(network.junctions.is_empty());
    }

    #[test]
    fn test_symmetric_tree_generation() {
        let network = BifurcationNetwork::<f64>::create_symmetric_tree(
            0.01,  // 1 cm root radius
            0.1,   // 10 cm root length
            2,     // 2 generations
            0.8,   // 80% length ratio
        );

        // 2 generations: 1 + 2 = 3 vessels minimum
        assert!(network.vessels.len() >= 3);
    }

    #[test]
    fn test_murray_validation() {
        let network = BifurcationNetwork::<f64>::create_symmetric_tree(
            0.01,
            0.1,
            2,
            0.8,
        );

        // Symmetric tree should satisfy Murray's Law
        let valid = network.validate_murrays_law(0.01);
        assert!(valid, "Symmetric tree should satisfy Murray's Law");
    }

    #[test]
    fn test_total_resistance() {
        let mut network = BifurcationNetwork::<f64>::new();
        network.add_vessel(VesselSegment::new(0, 0.003, 0.1, 0, 1));
        network.add_vessel(VesselSegment::new(1, 0.002, 0.05, 1, 2));

        let r = network.total_resistance();
        assert!(r > 0.0 && r.is_finite());
    }
}
