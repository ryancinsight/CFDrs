//! Flow analysis module for network systems.

use crate::channel::FlowRegime;
use nalgebra::RealField;
use std::collections::HashMap;

/// Comprehensive flow analysis for network systems
#[derive(Debug, Clone)]
pub struct FlowAnalysis<T: RealField + Copy> {
    /// Total flow rate through the network [m³/s]
    pub total_flow_rate: T,
    /// Flow rates through individual components [m³/s]
    pub component_flows: HashMap<String, T>,
    /// Average velocities in channels [m/s]
    pub velocities: HashMap<String, T>,
    /// Reynolds numbers for each channel
    pub reynolds_numbers: HashMap<String, T>,
    /// Flow regime classification
    pub flow_regimes: HashMap<String, FlowRegime>,
}

impl<T: RealField + Copy> FlowAnalysis<T> {
    /// Create a new flow analysis
    pub fn new() -> Self {
        Self {
            total_flow_rate: T::zero(),
            component_flows: HashMap::new(),
            velocities: HashMap::new(),
            reynolds_numbers: HashMap::new(),
            flow_regimes: HashMap::new(),
        }
    }

    /// Add flow data for a component
    pub fn add_component_flow(&mut self, id: String, flow_rate: T) {
        self.component_flows.insert(id, flow_rate);
        if flow_rate > T::zero() {
            self.total_flow_rate = self.total_flow_rate + flow_rate;
        }
    }

    /// Add velocity data for a component
    pub fn add_velocity(&mut self, id: String, velocity: T) {
        self.velocities.insert(id, velocity);
    }

    /// Add Reynolds number for a component
    pub fn add_reynolds_number(&mut self, id: String, reynolds: T) {
        self.reynolds_numbers.insert(id, reynolds);
    }

    /// Add flow regime classification for a component
    pub fn add_flow_regime(&mut self, id: String, regime: FlowRegime) {
        self.flow_regimes.insert(id, regime);
    }

    /// Get the average flow rate
    pub fn average_flow_rate(&self) -> T {
        if self.component_flows.is_empty() {
            T::zero()
        } else {
            let sum: T = self.component_flows.values().copied().sum();
            sum / T::from_usize(self.component_flows.len()).unwrap_or_else(T::one)
        }
    }

    /// Get the maximum flow rate
    pub fn max_flow_rate(&self) -> Option<T> {
        self.component_flows.values().max_by(|a, b| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        }).copied()
    }

    /// Get the minimum flow rate
    pub fn min_flow_rate(&self) -> Option<T> {
        self.component_flows.values().min_by(|a, b| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        }).copied()
    }
}

impl<T: RealField + Copy> Default for FlowAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}