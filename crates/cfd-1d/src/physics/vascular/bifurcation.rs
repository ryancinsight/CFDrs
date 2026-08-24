//! Vascular bifurcation network models
//!
//! This module provides 1D network models for arterial and venous bifurcations,
//! supporting hierarchical vascular trees with pressure-flow coupling.
//!
//! # Mathematical Foundation
//!
//! ## Theorem: 1D Bifurcation Energy Conservation
//!
//! **Theorem**: At a singular bifurcation junction connecting a parent vessel ($0$)
//! to $N$ daughter vessels ($i=1\dots N$), the steady-state 1D energy conservation
//! (Bernoulli equation integrated over the junction control volume) dictates:
//!
//! $$ p_0 + \frac{1}{2}\rho v_0^2 = p_i + \frac{1}{2}\rho v_i^2 + \Delta p_{\text{loss},i} \quad \forall i \in \{1 \dots N\} $$
//!
//! where the irreversible form loss $\Delta p_{\text{loss},i} = K_i \frac{1}{2}\rho v_0^2$.
//!
//! If kinetic energy recovery ($\frac{1}{2}\rho(v_0^2 - v_i^2)$) is assumed negligible
//! or explicitly captured by the empirical $K$-factor scaling, the junction simplifies
//! to a resistive node:
//! $$ p_i = p_0 - K_i \frac{1}{2}\rho v_0^2 $$
//!
//! Mass conservation is strictly simultaneously enforced: $Q_0 = \sum_{i=1}^N Q_i$.
//!
//! ## Network Resistance
//! Series: R_total = R₁ + R₂
//! Parallel: 1/R_total = 1/R₁ + 1/R₂
//!
//! # References
//! - Olufsen, M.S. (1999) "Structured tree outflow condition"
//! - Sherwin, S.J. et al. (2003) "One-dimensional modelling of a vascular network"

use super::murrays_law::MurraysLaw;
use crate::solver::analysis::resistance::{
    parallel_resistance as combine_parallel_resistances,
    series_resistance as combine_series_resistances,
};
use aequitas::systems::si::quantities::{
    Area, Compliance, DynamicViscosity, HydraulicInertance, HydraulicResistance, Length,
    MassDensity, Pressure, Velocity, VolumetricFlowRate,
};
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

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
pub struct VesselSegment<T: CfdScalar + Copy> {
    /// Segment identifier
    pub id: usize,
    /// Vessel radius \[m]
    pub radius: Length<T>,
    /// Vessel length \[m]
    pub length: Length<T>,
    /// Wall thickness \[m]
    pub wall_thickness: Length<T>,
    /// Young's modulus of wall \[Pa]
    pub youngs_modulus: Pressure<T>,
    /// Inlet node ID
    pub inlet_node: usize,
    /// Outlet node ID
    pub outlet_node: usize,
}

impl<T: CfdScalar + FloatElement + Copy> VesselSegment<T> {
    /// Create new vessel segment
    pub fn new(
        id: usize,
        radius: Length<T>,
        length: Length<T>,
        inlet_node: usize,
        outlet_node: usize,
    ) -> Self {
        let radius_base = radius.into_base();
        Self {
            id,
            radius,
            length,
            wall_thickness: Length::from_base(radius_base * <T as FloatElement>::from_f64(0.1)), // 10% of radius
            youngs_modulus: Pressure::from_base(<T as FloatElement>::from_f64(0.4e6)),
            // ~0.4 MPa for arteries
            inlet_node,
            outlet_node,
        }
    }

    /// Calculate Poiseuille resistance R = 8μL/(πR⁴)
    pub fn resistance(&self, viscosity: DynamicViscosity<T>) -> HydraulicResistance<T> {
        let pi = T::pi();
        let eight = <T as FloatElement>::from_f64(8.0);
        let viscosity = viscosity.into_base();
        let length = self.length.into_base();
        let radius = self.radius.into_base();
        HydraulicResistance::from_base(
            eight * viscosity * length / (pi * <T as FloatElement>::powi(radius, 4)),
        )
    }

    /// Calculate inertance L = ρL/(πR²)
    pub fn inertance(&self, density: MassDensity<T>) -> HydraulicInertance<T> {
        let pi = T::pi();
        let density = density.into_base();
        let length = self.length.into_base();
        let radius = self.radius.into_base();
        HydraulicInertance::from_base(density * length / (pi * radius * radius))
    }

    /// Calculate compliance C = 3πR³L/(2Eh)
    /// Based on thin-walled tube approximation
    pub fn compliance(&self) -> Compliance<T> {
        let pi = T::pi();
        let three = T::ONE + T::ONE + T::ONE;
        let two = T::ONE + T::ONE;
        let radius = self.radius.into_base();
        let length = self.length.into_base();
        let youngs_modulus = self.youngs_modulus.into_base();
        let wall_thickness = self.wall_thickness.into_base();
        Compliance::from_base(
            three * pi * <T as FloatElement>::powi(radius, 3) * length
                / (two * youngs_modulus * wall_thickness),
        )
    }

    /// Calculate wave speed c = √(Eh/(2ρR))
    pub fn wave_speed(&self, density: MassDensity<T>) -> Velocity<T> {
        let two = T::ONE + T::ONE;
        let density = density.into_base();
        let radius = self.radius.into_base();
        let wall_thickness = self.wall_thickness.into_base();
        let youngs_modulus = self.youngs_modulus.into_base();
        Velocity::from_base(<T as NumericElement>::sqrt(
            youngs_modulus * wall_thickness / (two * density * radius),
        ))
    }

    /// Calculate cross-sectional area
    pub fn area(&self) -> Area<T> {
        let pi = T::pi();
        let radius = self.radius.into_base();
        Area::from_base(pi * radius * radius)
    }

    /// Calculate diameter
    pub fn diameter(&self) -> Length<T> {
        let two = T::ONE + T::ONE;
        Length::from_base(two * self.radius.into_base())
    }
}

// ============================================================================
// Bifurcation Node
// ============================================================================

/// A bifurcation junction connecting vessel segments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bifurcation<T: CfdScalar + Copy> {
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
    /// Current pressure at junction \[Pa]
    pub pressure: Pressure<T>,
    /// Current flow rates [m³/s] - indexed by vessel ID
    pub flow_rates: Vec<VolumetricFlowRate<T>>,
}

impl<T: CfdScalar + FloatElement + Copy> Bifurcation<T> {
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
            k_factor: <T as FloatElement>::from_f64(0.1),
            parent_vessels: vec![parent_id],
            daughter_vessels: daughter_ids,
            pressure: Pressure::from_base(<T as FloatElement>::from_f64(13_332.0)),
            // ~100 mmHg
            flow_rates: vec![VolumetricFlowRate::from_base(T::ZERO); n_flows],
        }
    }

    /// Create confluence (venous) junction
    pub fn confluence(id: usize, parent_ids: Vec<usize>, daughter_id: usize) -> Self {
        let n_flows = parent_ids.len() + 1;
        Self {
            id,
            junction_type: JunctionType::Confluence,
            loss_model: JunctionLossModel::None,
            k_factor: <T as FloatElement>::from_f64(0.05),
            parent_vessels: parent_ids,
            daughter_vessels: vec![daughter_id],
            pressure: Pressure::from_base(<T as FloatElement>::from_f64(1_333.0)),
            // ~10 mmHg for veins
            flow_rates: vec![VolumetricFlowRate::from_base(T::ZERO); n_flows],
        }
    }

    /// Calculate mass conservation error
    pub fn mass_conservation_error(&self) -> T {
        if self.flow_rates.is_empty() {
            return T::ZERO;
        }

        let parent_flow: T = self
            .flow_rates
            .iter()
            .take(self.parent_vessels.len())
            .fold(T::ZERO, |acc, f| acc + f.into_base());
        let daughter_flow: T = self
            .flow_rates
            .iter()
            .skip(self.parent_vessels.len())
            .fold(T::ZERO, |acc, f| acc + f.into_base());

        if <T as NumericElement>::abs(parent_flow) < <T as FloatElement>::from_f64(1e-20) {
            return T::ZERO;
        }

        <T as NumericElement>::abs(parent_flow - daughter_flow)
            / <T as NumericElement>::abs(parent_flow)
    }

    /// Apply junction loss model to get pressure at daughters
    pub fn daughter_pressure(
        &self,
        parent_pressure: Pressure<T>,
        parent_velocity: Velocity<T>,
        density: MassDensity<T>,
    ) -> Pressure<T> {
        let parent_pressure_base = parent_pressure.into_base();
        let parent_velocity_base = parent_velocity.into_base();
        let density_base = density.into_base();
        match self.loss_model {
            JunctionLossModel::None => Pressure::from_base(parent_pressure_base),
            JunctionLossModel::KFactor => {
                let half = T::ONE / (T::ONE + T::ONE);
                Pressure::from_base(
                    parent_pressure_base
                        - self.k_factor
                            * half
                            * density_base
                            * parent_velocity_base
                            * parent_velocity_base,
                )
            }
            JunctionLossModel::EnergyPreserving => {
                // Energy preservation: p + ρv²/2 = const
                // More complex model would account for velocity changes
                Pressure::from_base(parent_pressure_base)
            }
        }
    }
}

// ============================================================================
// Bifurcation Network
// ============================================================================

/// A network of vessel segments connected by bifurcations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationNetwork<T: CfdScalar + Copy> {
    /// All vessel segments
    pub vessels: Vec<VesselSegment<T>>,
    /// All junction nodes
    pub junctions: Vec<Bifurcation<T>>,
    /// Inlet boundary pressure \[Pa]
    pub inlet_pressure: Pressure<T>,
    /// Outlet boundary resistance [Pa·s/m³]
    pub outlet_resistance: HydraulicResistance<T>,
    /// Blood density [kg/m³]
    pub density: MassDensity<T>,
    /// Blood viscosity [Pa·s]
    pub viscosity: DynamicViscosity<T>,
}

impl<T: CfdScalar + FloatElement + Copy> BifurcationNetwork<T> {
    /// Create empty network with default blood properties
    pub fn new() -> Self {
        Self {
            vessels: Vec::new(),
            junctions: Vec::new(),
            inlet_pressure: Pressure::from_base(<T as FloatElement>::from_f64(13_332.0)),
            // 100 mmHg
            outlet_resistance: HydraulicResistance::from_base(<T as FloatElement>::from_f64(1e9)),
            density: MassDensity::from_base(<T as FloatElement>::from_f64(
                cfd_core::physics::fluid::blood::constants::BLOOD_DENSITY,
            )),
            viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(
                cfd_core::physics::fluid::blood::constants::INFINITE_SHEAR_VISCOSITY,
            )),
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
    /// * `root_radius` - Root vessel radius \[m]
    /// * `root_length` - Root vessel length \[m]
    /// * `generations` - Number of bifurcation generations (1 = just root)
    /// * `length_ratio` - Length of daughter / length of parent
    pub fn create_symmetric_tree(
        root_radius: Length<T>,
        root_length: Length<T>,
        generations: usize,
        length_ratio: T,
    ) -> Self {
        let mut network = Self::new();
        let murray = MurraysLaw::<T>::new();

        // Recursive helper to build tree
        fn add_generation<T: CfdScalar + FloatElement + Copy>(
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
                Length::from_base(parent_radius),
                Length::from_base(parent_length),
                inlet_node,
                outlet_node,
            );
            network.vessels.push(vessel);

            if gen >= max_gen {
                return outlet_node;
            }

            // Create daughters
            let daughter_radius = murray
                .symmetric_daughter_diameter(parent_radius * (T::ONE + T::ONE) / (T::ONE + T::ONE));
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
            root_radius.into_base(),
            root_length.into_base(),
            length_ratio,
            1,
            generations,
            0,
        );

        network
    }

    /// Calculate total network resistance
    pub fn total_resistance(&self) -> HydraulicResistance<T> {
        if self.vessels.is_empty() {
            return HydraulicResistance::from_base(T::ZERO);
        }

        let mut parent_to_junction: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::with_capacity(self.junctions.len());
        let mut daughter_vessels = std::collections::HashSet::with_capacity(self.junctions.len());

        for (junction_index, junction) in self.junctions.iter().enumerate() {
            if let Some(&parent_vessel) = junction.parent_vessels.first() {
                parent_to_junction.insert(parent_vessel, junction_index);
            }

            for &daughter_vessel in &junction.daughter_vessels {
                daughter_vessels.insert(daughter_vessel);
            }
        }

        let mut root_vessels: Vec<usize> = self
            .vessels
            .iter()
            .filter(|vessel| !daughter_vessels.contains(&vessel.id))
            .map(|vessel| vessel.id)
            .collect();

        if root_vessels.is_empty() {
            return HydraulicResistance::from_base(T::ZERO);
        }

        root_vessels.sort_unstable();

        let mut memo = std::collections::HashMap::with_capacity(self.vessels.len());
        let mut root_resistances = Vec::with_capacity(root_vessels.len());

        for root_vessel in root_vessels {
            root_resistances.push(self.branch_equivalent_resistance(
                root_vessel,
                &parent_to_junction,
                &mut memo,
            ));
        }

        combine_parallel_resistances(
            root_resistances
                .into_iter()
                .map(HydraulicResistance::from_base),
        )
    }

    fn branch_equivalent_resistance(
        &self,
        vessel_id: usize,
        parent_to_junction: &std::collections::HashMap<usize, usize>,
        memo: &mut std::collections::HashMap<usize, T>,
    ) -> T {
        if let Some(&cached) = memo.get(&vessel_id) {
            return cached;
        }

        let Some(vessel) = self.vessels.get(vessel_id) else {
            return T::ZERO;
        };

        let series_resistance = vessel.resistance(self.viscosity).into_base();
        let equivalent = if let Some(&junction_index) = parent_to_junction.get(&vessel_id) {
            let junction = &self.junctions[junction_index];
            let mut downstream_resistances = Vec::with_capacity(junction.daughter_vessels.len());

            for &daughter_vessel in &junction.daughter_vessels {
                downstream_resistances.push(self.branch_equivalent_resistance(
                    daughter_vessel,
                    parent_to_junction,
                    memo,
                ));
            }

            combine_series_resistances([
                HydraulicResistance::from_base(series_resistance),
                combine_parallel_resistances(
                    downstream_resistances
                        .into_iter()
                        .map(HydraulicResistance::from_base),
                ),
            ])
        } else {
            HydraulicResistance::from_base(series_resistance)
        }
        .into_base();

        memo.insert(vessel_id, equivalent);
        equivalent
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

            if !murray.is_valid(
                parent.diameter().into_base(),
                d1.diameter().into_base(),
                d2.diameter().into_base(),
                tolerance,
            ) {
                return false;
            }
        }

        true
    }
}

impl<T: CfdScalar + FloatElement + Copy> Default for BifurcationNetwork<T> {
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

    #[test]
    fn test_vessel_segment_resistance() {
        let vessel =
            VesselSegment::<f64>::new(0, Length::from_base(0.003), Length::from_base(0.1), 0, 1);
        let r = vessel.resistance(DynamicViscosity::from_base(0.0035));

        // R = 8μL/(πR⁴) should be finite and positive
        assert!(r.into_base() > 0.0 && r.into_base().is_finite());
    }

    #[test]
    fn test_vessel_wave_speed() {
        let vessel =
            VesselSegment::<f64>::new(0, Length::from_base(0.003), Length::from_base(0.1), 0, 1);
        let c = vessel.wave_speed(MassDensity::from_base(1060.0));

        // Wave speed should be ~5-10 m/s for arteries
        assert!(
            c.into_base() > 1.0 && c.into_base() < 20.0,
            "Wave speed {} should be 5-10 m/s",
            c.into_base()
        );
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
            Length::from_base(0.01), // 1 cm root radius
            Length::from_base(0.1),  // 10 cm root length
            2,                       // 2 generations
            0.8,                     // 80% length ratio
        );

        // 2 generations: 1 + 2 = 3 vessels minimum
        assert!(network.vessels.len() >= 3);
    }

    #[test]
    fn test_murray_validation() {
        let network = BifurcationNetwork::<f64>::create_symmetric_tree(
            Length::from_base(0.01),
            Length::from_base(0.1),
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
        let root = VesselSegment::new(0, Length::from_base(0.003), Length::from_base(0.1), 0, 1);
        let daughter_a =
            VesselSegment::new(1, Length::from_base(0.002), Length::from_base(0.05), 1, 2);
        let daughter_b =
            VesselSegment::new(2, Length::from_base(0.002), Length::from_base(0.05), 1, 3);

        let root_resistance = root.resistance(network.viscosity).into_base();
        let daughter_a_resistance = daughter_a.resistance(network.viscosity).into_base();
        let daughter_b_resistance = daughter_b.resistance(network.viscosity).into_base();
        let expected =
            root_resistance + 1.0 / (1.0 / daughter_a_resistance + 1.0 / daughter_b_resistance);

        network.add_vessel(root);
        network.add_vessel(daughter_a);
        network.add_vessel(daughter_b);
        network.add_junction(Bifurcation::new(1, 0, vec![1, 2]));

        let r = network.total_resistance().into_base();
        assert!((r - expected).abs() < 1e-12, "got {r}, expected {expected}");
    }
}
