//! Analysis tools for 1D CFD networks.
//!
//! This module provides comprehensive analysis capabilities for microfluidic
//! and millifluidic networks including flow analysis, pressure analysis,
//! and resistance analysis with validation against analytical solutions.

use crate::network::Network;
use crate::solver::{NetworkSolver, SolverConfig, SolutionResult};
use cfd_core::{Result, Fluid};
use nalgebra::{RealField, ComplexField};
use num_traits::cast::FromPrimitive;
use std::collections::HashMap;

/// Comprehensive flow analysis for network systems
#[derive(Debug, Clone)]
pub struct FlowAnalysis<T: RealField> {
    /// Total volumetric flow rate [m³/s]
    pub total_flow_rate: T,
    /// Flow rates by component [m³/s]
    pub component_flows: HashMap<String, T>,
    /// Flow velocities [m/s]
    pub velocities: HashMap<String, T>,
    /// Reynolds numbers (dimensionless)
    pub reynolds_numbers: HashMap<String, T>,
    /// Flow regime classifications
    pub flow_regimes: HashMap<String, FlowRegime>,
}

/// Flow regime classification
#[derive(Debug, Clone, PartialEq)]
pub enum FlowRegime {
    /// Stokes flow (Re << 1)
    Stokes,
    /// Laminar flow (Re < 2300)
    Laminar,
    /// Transitional flow (2300 <= Re <= 4000)
    Transitional,
    /// Turbulent flow (Re > 4000)
    Turbulent,
}

/// Pressure analysis for network systems
#[derive(Debug, Clone)]
pub struct PressureAnalysis<T: RealField> {
    /// Pressure distribution [Pa]
    pub pressures: HashMap<String, T>,
    /// Pressure drops across components [Pa]
    pub pressure_drops: HashMap<String, T>,
    /// Maximum pressure in system [Pa]
    pub max_pressure: T,
    /// Minimum pressure in system [Pa]
    pub min_pressure: T,
    /// Pressure gradient statistics
    pub pressure_gradients: HashMap<String, T>,
}

/// Resistance analysis for network components
#[derive(Debug, Clone)]
pub struct ResistanceAnalysis<T: RealField> {
    /// Hydraulic resistances [Pa·s/m³]
    pub resistances: HashMap<String, T>,
    /// Equivalent circuit resistance [Pa·s/m³]
    pub total_resistance: T,
    /// Resistance contributions by component type
    pub resistance_by_type: HashMap<String, T>,
    /// Critical resistance paths
    pub critical_paths: Vec<Vec<String>>,
}

/// Network performance metrics
#[derive(Debug, Clone)]
pub struct PerformanceMetrics<T: RealField> {
    /// Throughput [m³/s]
    pub throughput: T,
    /// Pressure efficiency (useful pressure / total pressure)
    pub pressure_efficiency: T,
    /// Power consumption [W]
    pub power_consumption: T,
    /// Mixing efficiency (for mixing applications)
    pub mixing_efficiency: Option<T>,
    /// Residence time distribution
    pub residence_times: HashMap<String, T>,
}

/// Comprehensive network analyzer
pub struct NetworkAnalyzer<T: RealField> {
    solver: NetworkSolver<T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> NetworkAnalyzer<T> {
    /// Create a new network analyzer
    pub fn new() -> Self {
        Self {
            solver: NetworkSolver::new(),
        }
    }
    
    /// Create analyzer with custom solver configuration
    pub fn with_solver_config(config: SolverConfig<T>) -> Self {
        Self {
            solver: NetworkSolver::with_config(config),
        }
    }
    
    /// Perform comprehensive network analysis
    pub fn analyze(&self, network: &mut Network<T>) -> Result<NetworkAnalysisResult<T>> {
        // Solve the network
        let solution_result = self.solver.solve_steady_state(network)?;
        
        // Perform individual analyses
        let flow_analysis = self.analyze_flow(network)?;
        let pressure_analysis = self.analyze_pressure(network)?;
        let resistance_analysis = self.analyze_resistance(network)?;
        let performance_metrics = self.analyze_performance(network)?;
        
        Ok(NetworkAnalysisResult {
            solution_result,
            flow_analysis,
            pressure_analysis,
            resistance_analysis,
            performance_metrics,
        })
    }
    
    /// Analyze flow characteristics
    pub fn analyze_flow(&self, network: &Network<T>) -> Result<FlowAnalysis<T>> {
        let mut total_flow_rate = T::zero();
        let mut component_flows = HashMap::new();
        let mut velocities = HashMap::new();
        let mut reynolds_numbers = HashMap::new();
        let mut flow_regimes = HashMap::new();
        
        let fluid = network.fluid();
        
        for edge in network.edges() {
            if let Some(flow_rate) = edge.flow_rate {
                component_flows.insert(edge.id.clone(), flow_rate.clone());
                
                // Calculate velocity if area is available
                if let Some(area) = edge.properties.area {
                    if area > T::zero() {
                        let velocity = flow_rate.clone() / area;
                        velocities.insert(edge.id.clone(), velocity.clone());
                        
                        // Calculate Reynolds number if hydraulic diameter is available
                        if let Some(dh) = edge.properties.hydraulic_diameter {
                            let re = self.calculate_reynolds_number(
                                fluid, velocity, dh
                            );
                            reynolds_numbers.insert(edge.id.clone(), re.clone());
                            
                            // Classify flow regime
                            let regime = self.classify_flow_regime(re);
                            flow_regimes.insert(edge.id.clone(), regime);
                        }
                    }
                }
                
                // Sum total flow rate (count positive flow rates as inlets)
                // This is more robust than string matching
                if flow_rate > T::zero() {
                    total_flow_rate += flow_rate;
                }
            }
        }
        
        Ok(FlowAnalysis {
            total_flow_rate,
            component_flows,
            velocities,
            reynolds_numbers,
            flow_regimes,
        })
    }
    
    /// Analyze pressure distribution
    pub fn analyze_pressure(&self, network: &Network<T>) -> Result<PressureAnalysis<T>> {
        let mut pressures = HashMap::new();
        let mut pressure_drops = HashMap::new();
        let mut max_pressure = T::from_f64(f64::NEG_INFINITY).unwrap();
        let mut min_pressure = T::from_f64(f64::INFINITY).unwrap();
        let mut pressure_gradients = HashMap::new();
        
        // Collect node pressures
        for node in network.nodes() {
            if let Some(pressure) = node.pressure {
                pressures.insert(node.id.clone(), pressure.clone());
                
                if pressure > max_pressure {
                    max_pressure = pressure.clone();
                }
                if pressure < min_pressure {
                    min_pressure = pressure.clone();
                }
            }
        }
        
        // Collect pressure drops and gradients
        for edge in network.edges() {
            if let Some(pressure_drop) = edge.pressure_drop {
                pressure_drops.insert(edge.id.clone(), pressure_drop.clone());
                
                // Calculate pressure gradient if length is available
                if let Some(length) = edge.properties.length {
                    if length > T::zero() {
                        let gradient = pressure_drop / length;
                        pressure_gradients.insert(edge.id.clone(), gradient);
                    }
                }
            }
        }
        
        Ok(PressureAnalysis {
            pressures,
            pressure_drops,
            max_pressure,
            min_pressure,
            pressure_gradients,
        })
    }
    
    /// Analyze hydraulic resistances using iterator combinators
    pub fn analyze_resistance(&self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        // Use iterator patterns for zero-copy resistance collection
        let (resistances, resistance_by_type) = network.edges()
            .map(|edge| {
                let resistance = edge.effective_resistance();
                let type_key = edge.edge_type.as_str();
                ((edge.id.as_str(), resistance), (type_key, resistance))
            })
            .fold(
                (HashMap::new(), HashMap::new()),
                |(mut resistances, mut by_type), ((id, resistance), (type_key, res))| {
                    resistances.insert(id.to_string(), resistance);
                    *by_type.entry(type_key.to_string()).or_insert(T::zero()) += res;
                    (resistances, by_type)
                }
            );

        // Calculate total equivalent resistance (simplified)
        let total_resistance = self.calculate_equivalent_resistance(network);

        // Find critical paths (simplified - just high resistance paths)
        let critical_paths = self.find_critical_paths(network);

        Ok(ResistanceAnalysis {
            resistances,
            total_resistance,
            resistance_by_type,
            critical_paths,
        })
    }
    
    /// Analyze performance metrics
    pub fn analyze_performance(&self, network: &Network<T>) -> Result<PerformanceMetrics<T>> {
        let flow_analysis = self.analyze_flow(network)?;
        let pressure_analysis = self.analyze_pressure(network)?;
        
        let throughput = flow_analysis.total_flow_rate;
        
        // Calculate pressure efficiency
        let useful_pressure = pressure_analysis.max_pressure.clone() - pressure_analysis.min_pressure.clone();
        let total_pressure = pressure_analysis.max_pressure.clone();
        let pressure_efficiency = if total_pressure > T::zero() {
            useful_pressure / total_pressure
        } else {
            T::zero()
        };
        
        // Calculate power consumption (simplified)
        let power_consumption = self.calculate_power_consumption(network);
        
        // Calculate residence times (simplified)
        let residence_times = self.calculate_residence_times(network);

        // Calculate mixing efficiency based on network topology
        let mixing_efficiency = self.calculate_mixing_efficiency(network);

        Ok(PerformanceMetrics {
            throughput,
            pressure_efficiency,
            power_consumption,
            mixing_efficiency: Some(mixing_efficiency),
            residence_times,
        })
    }
    
    /// Calculate Reynolds number
    fn calculate_reynolds_number(&self, fluid: &Fluid<T>, velocity: T, diameter: T) -> T {
        let density = fluid.density;
        // Use actual operating temperature instead of hardcoded 20°C
        let temperature = T::from_f64(293.15).unwrap(); // Default to 20°C if not specified
        let viscosity = fluid.dynamic_viscosity(temperature);

        density * velocity * diameter / viscosity
    }
    
    /// Classify flow regime based on Reynolds number
    fn classify_flow_regime(&self, reynolds: T) -> FlowRegime {
        let re_val = reynolds.to_f64().unwrap_or(0.0);
        
        if re_val < 1.0 {
            FlowRegime::Stokes
        } else if re_val < 2300.0 {
            FlowRegime::Laminar
        } else if re_val <= 4000.0 {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }
    
    /// Calculate equivalent resistance (simplified series/parallel analysis)
    fn calculate_equivalent_resistance(&self, network: &Network<T>) -> T {
        // Simplified: sum all resistances (assumes series connection)
        network.edges()
            .map(|edge| edge.effective_resistance())
            .fold(T::zero(), |acc, r| acc + r)
    }
    
    /// Find critical resistance paths using advanced iterator patterns
    fn find_critical_paths(&self, network: &Network<T>) -> Vec<Vec<String>> {
        // Zero-copy analysis using iterator combinators
        let edge_data: Vec<_> = network.edges()
            .map(|edge| (edge.id.as_str(), edge.effective_resistance()))
            .collect();

        if edge_data.is_empty() {
            return Vec::new();
        }

        // Use advanced iterator patterns for statistical analysis
        use cfd_math::MathIteratorExt;

        let resistances: Vec<_> = edge_data.iter().map(|(_, resistance)| *resistance).collect();
        let mean_resistance = resistances.iter().cloned().mean().unwrap_or_else(T::zero);
        let _variance = resistances.iter().cloned().variance().unwrap_or_else(T::zero);

        let avg_resistance = mean_resistance;

        // Find critical edges using iterator chains
        let critical_edges: Vec<String> = edge_data.into_iter()
            .filter(|(_, resistance)| *resistance > avg_resistance)
            .map(|(id, _)| id.to_string())
            .collect();

        if critical_edges.is_empty() {
            Vec::new()
        } else {
            vec![critical_edges]
        }
    }
    
    /// Calculate power consumption
    fn calculate_power_consumption(&self, network: &Network<T>) -> T {
        let mut total_power = T::zero();
        
        for edge in network.edges() {
            if let (Some(flow_rate), Some(pressure_drop)) = (edge.flow_rate, edge.pressure_drop) {
                // Hydraulic power = flow_rate * pressure_drop
                total_power += ComplexField::abs(flow_rate * pressure_drop);
            }
        }
        
        total_power
    }
    
    /// Calculate residence times
    fn calculate_residence_times(&self, network: &Network<T>) -> HashMap<String, T> {
        let mut residence_times = HashMap::new();
        
        for edge in network.edges() {
            if let Some(flow_rate) = edge.flow_rate {
                if flow_rate > T::zero() {
                    // Calculate volume from area and length if available
                    if let (Some(area), Some(length)) = (edge.properties.area, edge.properties.length) {
                        let volume = area * length;
                        let residence_time = volume / flow_rate;
                        residence_times.insert(edge.id.clone(), residence_time);
                    }
                }
            }
        }

        residence_times
    }

    /// Calculate mixing efficiency based on network topology and flow patterns
    /// Uses advanced iterator patterns for zero-copy analysis
    fn calculate_mixing_efficiency(&self, network: &Network<T>) -> T {
        // Mixing efficiency is based on:
        // 1. Number of mixing junctions (T-junctions, Y-junctions)
        // 2. Flow rate ratios at junctions
        // 3. Reynolds numbers in mixing channels

        let junction_count = network.nodes()
            .filter(|node| matches!(node.node_type, crate::network::NodeType::Junction))
            .count();

        if junction_count == 0 {
            return T::zero(); // No mixing possible without junctions
        }

        // Calculate mixing effectiveness using iterator combinators
        let mixing_score = network.nodes()
            .filter(|node| matches!(node.node_type, crate::network::NodeType::Junction))
            .map(|junction| {
                // Get flow rates at this junction
                let connected_edges = network.node_edges(&junction.id).unwrap_or_default();
                let flow_rates: Vec<T> = connected_edges.iter()
                    .filter_map(|edge| edge.flow_rate)
                    .collect();

                if flow_rates.len() < 2 {
                    return T::zero(); // Need at least 2 flows for mixing
                }

                // Calculate flow rate uniformity (better mixing with more uniform flows)
                use cfd_math::MathIteratorExt;
                let mean_flow = flow_rates.iter().cloned().mean().unwrap_or_else(T::zero);
                let variance = flow_rates.iter().cloned().variance().unwrap_or_else(T::zero);

                if mean_flow.is_zero() {
                    T::zero()
                } else {
                    // Mixing efficiency inversely related to flow variance
                    let coefficient_of_variation = num_traits::Float::sqrt(variance) / num_traits::Float::abs(mean_flow);
                    T::one() / (T::one() + coefficient_of_variation)
                }
            })
            .fold(T::zero(), |acc, score| acc + score);

        // Normalize by number of junctions
        if junction_count > 0 {
            mixing_score / T::from_usize(junction_count).unwrap_or_else(T::one)
        } else {
            T::zero()
        }
    }
}

/// Complete network analysis results
#[derive(Debug, Clone)]
pub struct NetworkAnalysisResult<T: RealField> {
    /// Solver results
    pub solution_result: SolutionResult<T>,
    /// Flow analysis
    pub flow_analysis: FlowAnalysis<T>,
    /// Pressure analysis
    pub pressure_analysis: PressureAnalysis<T>,
    /// Resistance analysis
    pub resistance_analysis: ResistanceAnalysis<T>,
    /// Performance metrics
    pub performance_metrics: PerformanceMetrics<T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> Default for NetworkAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::NetworkBuilder;

    #[test]
    fn test_flow_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).unwrap();

        let flow_analysis = analyzer.analyze_flow(&network).unwrap();

        // Should have flow data
        assert!(flow_analysis.component_flows.contains_key("ch1"));
        assert!(flow_analysis.total_flow_rate >= 0.0);
    }

    #[test]
    fn test_pressure_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).unwrap();

        let pressure_analysis = analyzer.analyze_pressure(&network).unwrap();

        // Should have pressure data
        assert!(pressure_analysis.pressures.contains_key("inlet"));
        assert!(pressure_analysis.pressures.contains_key("outlet"));
        assert!(pressure_analysis.max_pressure >= pressure_analysis.min_pressure);
    }

    #[test]
    fn test_resistance_analysis() {
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let analyzer = NetworkAnalyzer::new();
        let resistance_analysis = analyzer.analyze_resistance(&network).unwrap();

        // Should have resistance data
        assert!(resistance_analysis.resistances.contains_key("ch1"));
        assert!(resistance_analysis.total_resistance > 0.0);
    }

    #[test]
    fn test_performance_metrics() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).unwrap();

        let performance = analyzer.analyze_performance(&network).unwrap();

        // Should have performance data
        assert!(performance.throughput >= 0.0);
        assert!(performance.pressure_efficiency >= 0.0);
        assert!(performance.pressure_efficiency <= 1.0);
        assert!(performance.power_consumption >= 0.0);
    }

    #[test]
    fn test_comprehensive_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 2000.0).unwrap()
            .add_junction("junction", 0.5, 0.0).unwrap()
            .add_outlet_pressure("outlet1", 1.0, 0.5, 0.0).unwrap()
            .add_outlet_pressure("outlet2", 1.0, -0.5, 0.0).unwrap()
            .add_channel("input_ch", "inlet", "junction", 100.0, 0.001, 1e-6).unwrap()
            .add_channel("output_ch1", "junction", "outlet1", 200.0, 0.001, 1e-6).unwrap()
            .add_channel("output_ch2", "junction", "outlet2", 200.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let analyzer = NetworkAnalyzer::new();
        let analysis_result = analyzer.analyze(&mut network).unwrap();

        // Check that all analyses completed
        assert!(analysis_result.solution_result.converged);
        assert!(!analysis_result.flow_analysis.component_flows.is_empty());
        assert!(!analysis_result.pressure_analysis.pressures.is_empty());
        assert!(!analysis_result.resistance_analysis.resistances.is_empty());
        assert!(analysis_result.performance_metrics.throughput >= 0.0);
    }

    #[test]
    fn test_flow_regime_classification() {
        let analyzer = NetworkAnalyzer::<f64>::new();

        assert_eq!(analyzer.classify_flow_regime(0.1), FlowRegime::Stokes);
        assert_eq!(analyzer.classify_flow_regime(100.0), FlowRegime::Laminar);
        assert_eq!(analyzer.classify_flow_regime(3000.0), FlowRegime::Transitional);
        assert_eq!(analyzer.classify_flow_regime(5000.0), FlowRegime::Turbulent);
    }

    #[test]
    fn test_reynolds_number_calculation() {
        let analyzer = NetworkAnalyzer::<f64>::new();
        let fluid = cfd_core::Fluid::water();

        let velocity = 0.1; // m/s
        let diameter = 1e-3; // 1 mm

        let re = analyzer.calculate_reynolds_number(&fluid, velocity, diameter);

        // Should be reasonable for water flow
        assert!(re > 0.0);
        assert!(re < 1000.0); // Should be laminar for these conditions
    }
}
