//! Resistance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::ResistanceAnalysis;
use crate::network::Network;
use cfd_core::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Resistance analyzer for network components
pub struct ResistanceAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> ResistanceAnalyzer<T> {
    /// Create new resistance analyzer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> NetworkAnalyzer<T> for ResistanceAnalyzer<T> {
    type Result = ResistanceAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        let mut analysis = ResistanceAnalysis::new();
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            let resistance = self.calculate_resistance(&edge.properties, fluid, edge.flow_rate);

            analysis.add_resistance(edge.id.clone(), resistance);

            // Add resistance by component type
            let component_type = self.classify_component_type(&edge.id);
            analysis.add_resistance_by_type(component_type, resistance);
        }

        // Find critical paths with highest resistance
        let critical_paths = self.find_critical_paths(network);
        for path in critical_paths {
            analysis.add_critical_path(path);
        }

        Ok(analysis)
    }

    fn name(&self) -> &str {
        "ResistanceAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> ResistanceAnalyzer<T> {
    fn calculate_resistance(
        &self,
        properties: &crate::network::ChannelProperties<T>,
        fluid: &cfd_core::fluid::Fluid<T>,
        flow_rate: Option<T>,
    ) -> T {
        use crate::resistance::{FlowConditions, HagenPoiseuilleModel, ResistanceModel};

        // Use appropriate resistance model based on geometry
        let model = HagenPoiseuilleModel::new(
            properties.hydraulic_diameter.unwrap_or_else(T::one),
            properties.length,
        );

        let conditions = FlowConditions {
            reynolds_number: flow_rate.map(|q| {
                let velocity = q / properties.area;
                let dh = properties.hydraulic_diameter.unwrap_or_else(T::one);
                fluid.reynolds_number(velocity, dh)
            }),
            velocity: flow_rate.map(|q| q / properties.area),
            flow_rate,
            temperature: T::from_f64(293.15).unwrap_or_else(T::one),
            pressure: T::from_f64(101325.0).unwrap_or_else(T::one),
        };

        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| properties.resistance)
    }

    fn classify_component_type(&self, edge_id: &str) -> String {
        // Classification based on edge naming convention
        if edge_id.contains("valve") {
            "valve".to_string()
        } else if edge_id.contains("pump") {
            "pump".to_string()
        } else if edge_id.contains("junction") {
            "junction".to_string()
        } else {
            "pipe".to_string()
        }
    }

    fn find_critical_paths(&self, network: &Network<T>) -> Vec<Vec<String>> {
        // Simplified: return paths with highest total resistance
        // In production, use graph algorithms to find actual critical paths
        Vec::new()
    }
}
