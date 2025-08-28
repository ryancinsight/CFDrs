//! Pressure analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::PressureAnalysis;
use crate::network::Network;
use cfd_core::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Pressure analyzer for network components
pub struct PressureAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> PressureAnalyzer<T> {
    /// Create new pressure analyzer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> NetworkAnalyzer<T> for PressureAnalyzer<T> {
    type Result = PressureAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<PressureAnalysis<T>> {
        let mut analysis = PressureAnalysis::new();

        // Collect node pressures
        let pressure_vec = network.pressures();
        for (idx, node) in network.nodes().enumerate() {
            if idx < pressure_vec.len() {
                let pressure = pressure_vec[idx];
                analysis.add_pressure(node.id.clone(), pressure);
            }
        }

        // Collect pressure drops and gradients
        for edge in network.edges_with_properties() {
            let (from_idx, to_idx) = edge.nodes;
            if from_idx < pressure_vec.len() && to_idx < pressure_vec.len() {
                let pressure_drop = pressure_vec[from_idx] - pressure_vec[to_idx];
                analysis.add_pressure_drop(edge.id.clone(), pressure_drop);

                // Calculate pressure gradient
                let length = edge.properties.length;
                if length > T::zero() {
                    let gradient = pressure_drop / length;
                    analysis.add_pressure_gradient(edge.id.clone(), gradient);
                }
            }
        }

        Ok(analysis)
    }

    fn name(&self) -> &str {
        "PressureAnalyzer"
    }
}
