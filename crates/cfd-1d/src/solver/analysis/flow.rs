//! Flow analysis module for network systems.

use super::blood_safety::{BloodShearLimits, HemolysisLimitViolation, ShearLimitViolation};
use crate::domain::channel::FlowRegime;
use crate::physics::hemolysis::{giersiepen_hi, taskin_hi};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use std::collections::HashMap;
use std::iter::Sum;

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
    /// Estimated wall shear rates for each channel [1/s]
    pub wall_shear_rates: HashMap<String, T>,
    /// Estimated wall shear stresses for each channel [Pa]
    pub wall_shear_stresses: HashMap<String, T>,
    /// Flow regime classification
    pub flow_regimes: HashMap<String, FlowRegime>,
}

impl<T: RealField + Copy + Sum> FlowAnalysis<T> {
    /// Create a new flow analysis
    #[must_use]
    pub fn new() -> Self {
        Self {
            total_flow_rate: T::zero(),
            component_flows: HashMap::new(),
            velocities: HashMap::new(),
            reynolds_numbers: HashMap::new(),
            wall_shear_rates: HashMap::new(),
            wall_shear_stresses: HashMap::new(),
            flow_regimes: HashMap::new(),
        }
    }

    /// Add flow data for a component
    pub fn add_component_flow(&mut self, id: String, flow_rate: T) {
        self.component_flows.insert(id, flow_rate);
        if flow_rate > T::zero() {
            self.total_flow_rate += flow_rate;
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

    /// Add wall shear rate for a component.
    pub fn add_wall_shear_rate(&mut self, id: String, shear_rate: T) {
        self.wall_shear_rates.insert(id, shear_rate);
    }

    /// Add wall shear stress for a component.
    pub fn add_wall_shear_stress(&mut self, id: String, shear_stress: T) {
        self.wall_shear_stresses.insert(id, shear_stress);
    }

    /// Add flow regime classification for a component
    pub fn add_flow_regime(&mut self, id: String, regime: FlowRegime) {
        self.flow_regimes.insert(id, regime);
    }

    /// Set total flow rate
    pub fn set_total_flow(&mut self, total_flow: T) {
        self.total_flow_rate = total_flow;
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
        self.component_flows
            .values()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
    }

    /// Get the minimum flow rate
    pub fn min_flow_rate(&self) -> Option<T> {
        self.component_flows
            .values()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
    }

    /// Flag components that exceed configured blood shear limits.
    #[must_use]
    pub fn flag_fda_shear_limit_violations(
        &self,
        limits: &BloodShearLimits<T>,
    ) -> Vec<ShearLimitViolation<T>> {
        let mut violations = Vec::new();

        for (component_id, wall_shear_stress) in &self.wall_shear_stresses {
            let exceeds_stress = *wall_shear_stress > limits.max_wall_shear_stress_pa;
            let shear_rate = self.wall_shear_rates.get(component_id).copied();
            let exceeds_rate = if let (Some(rate), Some(rate_limit)) =
                (shear_rate, limits.max_wall_shear_rate_per_s)
            {
                rate > rate_limit
            } else {
                false
            };

            if exceeds_stress || exceeds_rate {
                let ratio = *wall_shear_stress / limits.max_wall_shear_stress_pa;
                violations.push(ShearLimitViolation {
                    component_id: component_id.clone(),
                    wall_shear_stress_pa: *wall_shear_stress,
                    stress_limit_pa: limits.max_wall_shear_stress_pa,
                    stress_exceedance_ratio: ratio,
                    wall_shear_rate_per_s: shear_rate,
                    shear_rate_limit_per_s: limits.max_wall_shear_rate_per_s,
                });
            }
        }

        violations
    }
}

impl<T: RealField + Copy + Sum + FromPrimitive + ToPrimitive> FlowAnalysis<T> {
    /// Flag components whose time-integrated hemolysis exceeds configured limits.
    ///
    /// Residence times are provided externally because the reduced-order flow
    /// analysis stores local stress and velocity, but not a unique transit-time
    /// model for every downstream biological payload.
    #[must_use]
    pub fn flag_hemolysis_limit_violations(
        &self,
        limits: &BloodShearLimits<T>,
        residence_times_s: &HashMap<String, T>,
    ) -> Vec<HemolysisLimitViolation<T>> {
        if limits.max_giersiepen_hi.is_none() && limits.max_taskin_hi.is_none() {
            return Vec::new();
        }

        let mut violations = Vec::new();

        for (component_id, wall_shear_stress) in &self.wall_shear_stresses {
            let Some(exposure_time_s) = residence_times_s.get(component_id).copied() else {
                continue;
            };

            let shear_pa = wall_shear_stress.to_f64().unwrap_or(0.0);
            let duration_s = exposure_time_s.to_f64().unwrap_or(0.0);

            let giersiepen_value = limits
                .max_giersiepen_hi
                .and_then(|_| T::from_f64(giersiepen_hi(shear_pa, duration_s)));
            let taskin_value = limits
                .max_taskin_hi
                .and_then(|_| T::from_f64(taskin_hi(shear_pa, duration_s)));

            let giersiepen_ratio = match (giersiepen_value, limits.max_giersiepen_hi) {
                (Some(value), Some(limit)) if limit > T::zero() => Some(value / limit),
                _ => None,
            };
            let taskin_ratio = match (taskin_value, limits.max_taskin_hi) {
                (Some(value), Some(limit)) if limit > T::zero() => Some(value / limit),
                _ => None,
            };

            let exceeds_giersiepen = matches!(
                (giersiepen_value, limits.max_giersiepen_hi),
                (Some(value), Some(limit)) if value > limit
            );
            let exceeds_taskin = matches!(
                (taskin_value, limits.max_taskin_hi),
                (Some(value), Some(limit)) if value > limit
            );

            if exceeds_giersiepen || exceeds_taskin {
                violations.push(HemolysisLimitViolation {
                    component_id: component_id.clone(),
                    wall_shear_stress_pa: *wall_shear_stress,
                    exposure_time_s,
                    giersiepen_hi: giersiepen_value,
                    giersiepen_limit: limits.max_giersiepen_hi,
                    giersiepen_exceedance_ratio: giersiepen_ratio,
                    taskin_hi: taskin_value,
                    taskin_limit: limits.max_taskin_hi,
                    taskin_exceedance_ratio: taskin_ratio,
                });
            }
        }

        violations
    }
}

impl<T: RealField + Copy + Sum> Default for FlowAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::analysis::BloodShearLimits;

    #[test]
    fn flags_component_when_stress_exceeds_limit() {
        let mut analysis = FlowAnalysis::<f64>::new();
        analysis.add_wall_shear_stress("edge_1".to_string(), 180.0);
        analysis.add_wall_shear_rate("edge_1".to_string(), 40_000.0);

        let limits = BloodShearLimits {
            max_wall_shear_stress_pa: 150.0,
            max_wall_shear_rate_per_s: None,
            max_giersiepen_hi: None,
            max_taskin_hi: None,
        };

        let violations = analysis.flag_fda_shear_limit_violations(&limits);
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].component_id, "edge_1");
        assert!(violations[0].stress_exceedance_ratio > 1.0);
    }

    #[test]
    fn flags_component_when_rate_exceeds_optional_limit() {
        let mut analysis = FlowAnalysis::<f64>::new();
        analysis.add_wall_shear_stress("edge_2".to_string(), 80.0);
        analysis.add_wall_shear_rate("edge_2".to_string(), 50_000.0);

        let limits = BloodShearLimits {
            max_wall_shear_stress_pa: 150.0,
            max_wall_shear_rate_per_s: Some(40_000.0),
            max_giersiepen_hi: None,
            max_taskin_hi: None,
        };

        let violations = analysis.flag_fda_shear_limit_violations(&limits);
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].component_id, "edge_2");
    }

    #[test]
    fn flags_component_when_hemolysis_exceeds_limit() {
        let mut analysis = FlowAnalysis::<f64>::new();
        analysis.add_wall_shear_stress("edge_3".to_string(), 180.0);

        let limits = BloodShearLimits::fda_conservative_whole_blood()
            .with_hemolysis_limits(Some(1e-3), Some(5e-3));
        let residence_times = HashMap::from([("edge_3".to_string(), 0.5)]);

        let violations = analysis.flag_hemolysis_limit_violations(&limits, &residence_times);

        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].component_id, "edge_3");
        assert!(violations[0]
            .giersiepen_exceedance_ratio
            .is_some_and(|ratio| ratio > 1.0));
        assert!(violations[0]
            .taskin_exceedance_ratio
            .is_some_and(|ratio| ratio > 1.0));
    }
}
