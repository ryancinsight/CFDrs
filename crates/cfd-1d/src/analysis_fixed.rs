//! Fixed analysis module - temporary file for fixing compilation errors

use crate::network::{Network, EdgeProperties};
use nalgebra::{RealField, ComplexField};
use std::collections::HashMap;
use num_traits::FromPrimitive;

// Fixed calculate_power_consumption function
pub fn calculate_power_consumption<T: RealField + FromPrimitive + Copy>(
    network: &Network<T>
) -> T {
    let mut total_power = T::zero();
    let pressure_vec = network.pressures();
    
    for edge in network.edges_with_properties() {
        if let Some(flow_rate) = edge.flow_rate {
            let (from_idx, to_idx) = edge.nodes;
            if from_idx < pressure_vec.len() && to_idx < pressure_vec.len() {
                let pressure_drop = pressure_vec[from_idx] - pressure_vec[to_idx];
                // Hydraulic power = flow_rate * pressure_drop
                total_power += ComplexField::abs(flow_rate * pressure_drop);
            }
        }
    }
    
    total_power
}

// Fixed calculate_residence_times function
pub fn calculate_residence_times<T: RealField + FromPrimitive + Copy>(
    network: &Network<T>
) -> HashMap<String, T> {
    let mut residence_times = HashMap::new();
    
    for edge in network.edges_with_properties() {
        if let Some(flow_rate) = edge.flow_rate {
            if flow_rate > T::zero() {
                // Calculate volume from area and length (not Options)
                let area = edge.properties.area;
                let length = edge.properties.length;
                if area > T::zero() && length > T::zero() {
                    let volume = area * length;
                    let residence_time = volume / flow_rate;
                    residence_times.insert(edge.id.clone(), residence_time);
                }
            }
        }
    }

    residence_times
}