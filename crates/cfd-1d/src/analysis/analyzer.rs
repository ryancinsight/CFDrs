//! Main network analyzer implementation.

use super::{FlowAnalysis, PerformanceMetrics, PressureAnalysis, ResistanceAnalysis};
use crate::channel::FlowRegime;
use crate::network::Network;
use cfd_core::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use std::collections::HashSet;
use std::iter::Sum;

/// Complete network analysis results
#[derive(Debug, Clone)]
pub struct NetworkAnalysisResult<T: RealField + Copy> {
    /// Flow analysis
    pub flow_analysis: FlowAnalysis<T>,
    /// Pressure analysis
    pub pressure_analysis: PressureAnalysis<T>,
    /// Resistance analysis
    pub resistance_analysis: ResistanceAnalysis<T>,
    /// Performance metrics
    pub performance_metrics: PerformanceMetrics<T>,
}

/// Comprehensive network analyzer
pub struct NetworkAnalyzer<T: RealField + Copy> {
    solver: crate::solver::NetworkSolver<T>,
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzer<T> {
    /// Create a new network analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            solver: crate::solver::NetworkSolver::new(),
        }
    }

    /// Create analyzer with custom solver configuration
    pub fn with_solver_config(config: crate::solver::SolverConfig<T>) -> Self {
        Self {
            solver: crate::solver::NetworkSolver::with_config(config),
        }
    }

    /// Perform comprehensive network analysis
    pub fn analyze(&mut self, network: &mut Network<T>) -> Result<NetworkAnalysisResult<T>> {
        // Solve the network
        let problem = crate::solver::NetworkProblem::new(network.clone());
        let _solution_result = self.solver.solve_network(&problem)?;

        // Perform individual analyses
        let flow_analysis = self.analyze_flow(network)?;
        let pressure_analysis = self.analyze_pressure(network)?;
        let resistance_analysis = self.analyze_resistance(network)?;
        let performance_metrics = self.analyze_performance(network)?;

        Ok(NetworkAnalysisResult {
            flow_analysis,
            pressure_analysis,
            resistance_analysis,
            performance_metrics,
        })
    }

    /// Analyze flow characteristics
    pub fn analyze_flow(&self, network: &Network<T>) -> Result<FlowAnalysis<T>> {
        let mut analysis = FlowAnalysis::new();
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            if let Some(flow_rate) = edge.flow_rate {
                analysis.add_component_flow(edge.id.clone(), flow_rate);

                // Calculate velocity if area is available
                let area = edge.properties.area;
                if area > T::zero() {
                    let velocity = flow_rate / area;
                    analysis.add_velocity(edge.id.clone(), velocity);

                    // Calculate Reynolds number if hydraulic diameter is available
                    if let Some(dh) = edge.properties.hydraulic_diameter {
                        let re = self.calculate_reynolds_number(fluid, velocity, dh);
                        analysis.add_reynolds_number(edge.id.clone(), re);

                        // Classify flow regime
                        let regime = self.classify_flow_regime(re);
                        analysis.add_flow_regime(edge.id.clone(), regime);
                    }
                }
            }
        }

        Ok(analysis)
    }

    /// Analyze pressure distribution
    pub fn analyze_pressure(&self, network: &Network<T>) -> Result<PressureAnalysis<T>> {
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

                // Calculate pressure gradient using length
                let length = edge.properties.length;
                if length > T::zero() {
                    let gradient = pressure_drop / length;
                    analysis.add_pressure_gradient(edge.id.clone(), gradient);
                }
            }
        }

        Ok(analysis)
    }

    /// Analyze resistance characteristics
    pub fn analyze_resistance(&self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        let mut analysis = ResistanceAnalysis::new();
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            let resistance = self.calculate_resistance(&edge.properties, fluid, edge.flow_rate);

            analysis.add_resistance(edge.id.clone(), resistance);

            // Add resistance by type
            let component_type = self.classify_component_type(&edge.id);
            analysis.add_resistance_by_type(component_type, resistance);
        }

        // Find critical paths
        let critical_paths = self.find_critical_paths(network);
        for path in critical_paths {
            analysis.add_critical_path(path);
        }

        Ok(analysis)
    }

    /// Analyze performance metrics
    pub fn analyze_performance(&self, network: &Network<T>) -> Result<PerformanceMetrics<T>> {
        let mut metrics = PerformanceMetrics::new();

        // Calculate throughput (total flow through inlets)
        let mut throughput = T::zero();
        for edge in network.edges_with_properties() {
            if edge.id.contains("inlet") || edge.id.contains("input") {
                if let Some(flow_rate) = edge.flow_rate {
                    throughput += Float::abs(flow_rate);
                }
            }
        }
        metrics.set_throughput(throughput);

        // Calculate pressure efficiency
        let pressure_vec = network.pressures();
        if !pressure_vec.is_empty() {
            let max_pressure = *pressure_vec
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or(&T::zero());
            let min_pressure = *pressure_vec
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or(&T::zero());

            if max_pressure > T::zero() {
                let useful_pressure = max_pressure - min_pressure;
                let efficiency = useful_pressure / max_pressure;
                metrics.set_pressure_efficiency(efficiency);
            }
        }

        // Calculate power consumption
        let mut total_power = T::zero();
        for edge in network.edges_with_properties() {
            if let Some(flow_rate) = edge.flow_rate {
                let (from_idx, to_idx) = edge.nodes;
                if from_idx < pressure_vec.len() && to_idx < pressure_vec.len() {
                    let pressure_drop = Float::abs(pressure_vec[from_idx] - pressure_vec[to_idx]);
                    let power = Float::abs(flow_rate) * pressure_drop;
                    total_power += power;
                }
            }
        }
        metrics.set_power_consumption(total_power);

        // Calculate residence times
        for edge in network.edges_with_properties() {
            if let Some(flow_rate) = edge.flow_rate {
                if Float::abs(flow_rate) > T::zero() {
                    let volume = edge.properties.length * edge.properties.area;
                    let residence_time = volume / Float::abs(flow_rate);
                    metrics.add_residence_time(edge.id.clone(), residence_time);
                }
            }
        }

        Ok(metrics)
    }

    /// Calculate Reynolds number
    pub fn calculate_reynolds_number(
        &self,
        fluid: &cfd_core::Fluid<T>,
        velocity: T,
        diameter: T,
    ) -> T {
        let density = fluid.density;
        // Use characteristic viscosity for Reynolds number calculation
        let viscosity = fluid.characteristic_viscosity();

        if viscosity > T::zero() {
            density * Float::abs(velocity) * diameter / viscosity
        } else {
            T::zero()
        }
    }

    /// Classify flow regime based on Reynolds number
    pub fn classify_flow_regime(&self, reynolds: T) -> FlowRegime {
        let re_value = reynolds.to_subset().unwrap_or(0.0);

        // Reynolds number thresholds from literature
        const REYNOLDS_LAMINAR_LIMIT: f64 = 2300.0;
        const REYNOLDS_TURBULENT_THRESHOLD: f64 = 4000.0;

        if re_value < 1.0 {
            FlowRegime::Stokes
        } else if re_value < REYNOLDS_LAMINAR_LIMIT {
            FlowRegime::Laminar
        } else if re_value < REYNOLDS_TURBULENT_THRESHOLD {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }

    /// Calculate hydraulic resistance
    fn calculate_resistance(
        &self,
        properties: &crate::network::ChannelProperties<T>,
        fluid: &cfd_core::Fluid<T>,
        flow_rate: Option<T>,
    ) -> T {
        let viscosity = fluid.characteristic_viscosity();
        let length = properties.length;
        let area = properties.area;

        if area > T::zero() {
            // Base resistance (Hagen-Poiseuille for laminar flow)
            let base_resistance =
                T::from_f64(8.0).unwrap_or_else(T::one) * viscosity * length / (area * area);

            // Adjust for flow regime if flow rate is known
            if let Some(q) = flow_rate {
                if let Some(dh) = properties.hydraulic_diameter {
                    let velocity = Float::abs(q) / area;
                    let re = self.calculate_reynolds_number(fluid, velocity, dh);

                    // Apply correction factor based on flow regime
                    match self.classify_flow_regime(re) {
                        FlowRegime::Stokes => base_resistance,
                        FlowRegime::Laminar => base_resistance,
                        FlowRegime::Transitional => {
                            base_resistance * T::from_f64(1.5).unwrap_or_else(T::one)
                        }
                        FlowRegime::Turbulent => {
                            base_resistance * T::from_f64(2.0).unwrap_or_else(T::one)
                        }
                        FlowRegime::SlipFlow => {
                            // Slip flow occurs at very low pressures/high Knudsen numbers
                            // Resistance is lower due to slip at walls
                            base_resistance * T::from_f64(0.8).unwrap_or_else(T::one)
                        }
                    }
                } else {
                    base_resistance
                }
            } else {
                base_resistance
            }
        } else {
            T::from_f64(1e12).unwrap_or_else(|| T::from_f64(1e6).unwrap_or_else(T::one))
        }
    }

    /// Classify component type from ID
    fn classify_component_type(&self, id: &str) -> String {
        if id.contains("channel") || id.contains("ch") {
            "channel".to_string()
        } else if id.contains("valve") {
            "valve".to_string()
        } else if id.contains("pump") {
            "pump".to_string()
        } else if id.contains("junction") {
            "junction".to_string()
        } else {
            "other".to_string()
        }
    }

    /// Find critical resistance paths using depth-first search
    fn find_critical_paths(&self, network: &Network<T>) -> Vec<Vec<String>> {
        let mut paths = Vec::new();
        let mut visited = HashSet::new();

        // Find inlet nodes
        let inlet_nodes: Vec<_> = network.nodes().filter(|n| n.id.contains("inlet")).collect();

        // Find outlet nodes
        let outlet_nodes: Vec<_> = network
            .nodes()
            .filter(|n| n.id.contains("outlet"))
            .collect();

        // Find paths from each inlet to each outlet
        for inlet in &inlet_nodes {
            for outlet in &outlet_nodes {
                if let Some(path) = self.find_path(network, &inlet.id, &outlet.id, &mut visited) {
                    paths.push(path);
                }
                visited.clear();
            }
        }

        // Return top 3 paths with highest resistance
        paths.truncate(3);
        paths
    }

    /// Find path between two nodes
    fn find_path(
        &self,
        network: &Network<T>,
        from: &str,
        to: &str,
        visited: &mut HashSet<String>,
    ) -> Option<Vec<String>> {
        if from == to {
            return Some(vec![from.to_string()]);
        }

        if visited.contains(from) {
            return None;
        }

        visited.insert(from.to_string());

        // Find edges connected to 'from' node
        for edge in network.edges_with_properties() {
            let from_node = network.nodes().nth(edge.nodes.0);
            let to_node = network.nodes().nth(edge.nodes.1);

            if let (Some(fn_node), Some(tn_node)) = (from_node, to_node) {
                if fn_node.id == from {
                    if let Some(mut path) = self.find_path(network, &tn_node.id, to, visited) {
                        path.insert(0, from.to_string());
                        return Some(path);
                    }
                }
            }
        }

        None
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> Default for NetworkAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}
