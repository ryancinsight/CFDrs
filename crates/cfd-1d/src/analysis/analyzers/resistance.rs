//! Resistance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::error::ResistanceCalculationError;
use crate::analysis::ResistanceAnalysis;
use crate::network::Network;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use std::iter::Sum;

/// Resistance analyzer for network components
pub struct ResistanceAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for ResistanceAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> ResistanceAnalyzer<T> {
    /// Create new resistance analyzer
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzer<T>
    for ResistanceAnalyzer<T>
{
    type Result = ResistanceAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        let mut analysis = ResistanceAnalysis::new();
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            let flow_rate = if edge.flow_rate == T::zero() {
                None
            } else {
                Some(edge.flow_rate)
            };

            // Calculate resistance with proper error handling
            let resistance = self
                .calculate_resistance(edge.properties, fluid, flow_rate)
                .map_err(|e| {
                    // Enhance error context before propagating
                    cfd_core::error::Error::InvalidInput(format!(
                        "Failed to analyze resistance for edge '{}': {}",
                        edge.id, e
                    ))
                })?;

            analysis.add_resistance(edge.id.clone(), resistance);

            // Add resistance by component type (now type-safe)
            let component_type = edge.properties.component_type;
            analysis.add_resistance_by_type(component_type.as_str().to_string(), resistance);
        }

        // Note: Critical path analysis removed as it was unimplemented
        // This feature should be added as a separate, properly implemented method
        // when the algorithm is ready

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "ResistanceAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> ResistanceAnalyzer<T> {
    fn calculate_resistance(
        &self,
        properties: &crate::network::EdgeProperties<T>,
        fluid: &cfd_core::fluid::Fluid<T>,
        flow_rate: Option<T>,
    ) -> std::result::Result<T, ResistanceCalculationError> {
        use crate::resistance::{FlowConditions, HagenPoiseuilleModel, ResistanceModel};

        // Require hydraulic diameter - no silent fallbacks
        let hydraulic_diameter = properties
            .hydraulic_diameter
            .ok_or(ResistanceCalculationError::MissingHydraulicDiameter)?;

        // Create resistance model with validated parameters
        let model = HagenPoiseuilleModel::new(hydraulic_diameter, properties.length);

        // Build flow conditions with proper constants
        let conditions = FlowConditions {
            reynolds_number: flow_rate.map(|q| {
                let velocity = q / properties.area;
                fluid.density * velocity * hydraulic_diameter / fluid.viscosity
            }),
            velocity: flow_rate.map(|q| q / properties.area),
            flow_rate,
            // Use expect with clear messages for constants that must succeed
            temperature: T::from_f64(293.15)
                .expect("Standard temperature (293.15K) must be representable in type T"),
            pressure: T::from_f64(101_325.0)
                .expect("Standard pressure (101_325.0 Pa) must be representable in type T"),
        };

        // Calculate resistance and propagate any errors
        model
            .calculate_resistance(fluid, &conditions)
            .map_err(|e| ResistanceCalculationError::ModelError(e.to_string()))
    }
}
