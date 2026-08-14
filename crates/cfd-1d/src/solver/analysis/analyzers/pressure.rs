//! Pressure analysis for network components

use super::traits::NetworkAnalyzer;
use crate::domain::network::Network;
use crate::solver::analysis::PressureAnalysis;
use aequitas::systems::si::quantities::{Pressure, PressureGradient};
use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use cfd_core::error::Result;
use cfd_core::CfdScalar;
use std::iter::Sum;

/// Pressure analyzer for network components
pub struct PressureAnalyzer<T: CfdScalar + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: CfdScalar + Copy> Default for PressureAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: CfdScalar + Copy> PressureAnalyzer<T> {
    /// Create new pressure analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: CfdScalar + Copy + SafeFromF64 + SafeFromUsize + Sum> NetworkAnalyzer<T>
    for PressureAnalyzer<T>
{
    type Result = PressureAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<PressureAnalysis<T>> {
        let mut analysis = PressureAnalysis::new();

        // Collect node pressures

        let pressures = network.pressures();
        for (idx, node) in network.nodes().enumerate() {
            if let Some(&pressure) = pressures.get(idx) {
                analysis.add_pressure(node.id.clone(), Pressure::from_base(pressure.into_base()));
            }
        }

        // Collect pressure drops and gradients
        for edge in network.edges_with_properties() {
            let (from_idx, to_idx) = edge.nodes;
            if let (Some(&p_from), Some(&p_to)) = (
                pressures.get(from_idx.index()),
                pressures.get(to_idx.index()),
            ) {
                let pressure_drop = p_from.into_base() - p_to.into_base();
                analysis.add_pressure_drop(edge.id.clone(), Pressure::from_base(pressure_drop));

                // Calculate pressure gradient
                let length = edge.properties.length.into_base();
                if length > T::ZERO {
                    let gradient = pressure_drop / length;
                    analysis.add_pressure_gradient(
                        edge.id.clone(),
                        PressureGradient::from_base(gradient),
                    );
                }
            }
        }

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "PressureAnalyzer"
    }
}
