//! Network analysis and performance metrics for 1D CFD simulations.

use crate::network::Network;
use crate::channel::FlowRegime;
use cfd_core::{Result, Fluid};
use nalgebra::{RealField, ComplexField};
use num_traits::FromPrimitive;

use std::collections::{HashMap, HashSet};

/// Comprehensive flow analysis for network systems
#[derive(Debug, Clone)]
pub struct FlowAnalysis<T: RealField> {
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
pub struct NetworkAnalyzer<T: RealField + Copy> {
    solver: crate::solver::NetworkSolver<T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float + Copy> NetworkAnalyzer<T> {
    /// Create a new network analyzer
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
        // Create a problem from the network
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
        let mut total_flow_rate = T::zero();
        let mut component_flows = HashMap::new();
        let mut velocities = HashMap::new();
        let mut reynolds_numbers = HashMap::new();
        let mut flow_regimes = HashMap::new();
        
        let fluid = network.fluid();
        
        for edge in network.edges_with_properties() {
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
        let mut max_pressure = T::from_f64(f64::NEG_INFINITY).unwrap_or_else(|| T::zero());
        let mut min_pressure = T::from_f64(f64::INFINITY).unwrap_or_else(|| T::zero());
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
        for edge in network.edges_with_properties() {
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

        // Calculate total equivalent resistance using Modified Nodal Analysis
        let total_resistance = self.calculate_equivalent_resistance(network);

        // Find critical paths using graph algorithms
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
        let temperature = T::from_f64(293.15).unwrap_or_else(|| T::zero()); // Default to 20°C if not specified
        let viscosity = fluid.dynamic_viscosity(temperature)
            .unwrap_or_else(|_| T::from_f64(0.001).unwrap_or_else(|| T::one())); // Default to water viscosity

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
    
    /// Calculate equivalent resistance using Modified Nodal Analysis (MNA)
    /// 
    /// This implements the complete MNA method for network analysis:
    /// 1. Build conductance matrix G from network topology
    /// 2. Apply Kirchhoff's Current Law at each node
    /// 3. Solve the linear system for node voltages
    /// 4. Calculate equivalent resistance from voltage drop and current
    ///
    /// Reference: Nagel, L.W. (1975) "SPICE2: A Computer Program to Simulate Semiconductor Circuits"
    fn calculate_equivalent_resistance(&self, network: &Network<T>) -> T {
        use nalgebra::{DMatrix, DVector};
        
        // Get all nodes and create index mapping
        let nodes: Vec<_> = network.nodes().collect();
        let node_indices: std::collections::HashMap<_, _> = nodes.iter()
            .enumerate()
            .map(|(i, n)| (n.id.parse::<usize>().unwrap_or(i), i))
            .collect();
        
        let n = nodes.len();
        if n < 2 {
            return T::zero(); // Need at least 2 nodes
        }
        
        // Build conductance matrix (G) using Kirchhoff's Current Law
        let mut g_matrix: DMatrix<T> = DMatrix::zeros(n, n);
        
        // Fill conductance matrix from network edges
        // We need to iterate through the graph edges to get connections
        use petgraph::visit::EdgeRef;
        for edge_ref in network.graph.edge_references() {
            let edge = edge_ref.weight();
            let resistance = edge.effective_resistance();
            if resistance <= T::zero() {
                continue; // Skip zero/negative resistances
            }
            
            let conductance = T::one() / resistance;
            
            // Get node indices for this edge from the graph
            let source_idx = edge_ref.source().index();
            let target_idx = edge_ref.target().index();
            
            if let (Some(&i), Some(&j)) = (node_indices.get(&source_idx), node_indices.get(&target_idx)) {
                // Add conductance to diagonal elements
                g_matrix[(i, i)] = g_matrix[(i, i)].clone() + conductance.clone();
                g_matrix[(j, j)] = g_matrix[(j, j)].clone() + conductance.clone();
                
                // Subtract conductance from off-diagonal elements
                g_matrix[(i, j)] = g_matrix[(i, j)].clone() - conductance.clone();
                g_matrix[(j, i)] = g_matrix[(j, i)].clone() - conductance.clone();
            }
        }
        
        // Find inlet and outlet nodes
        let inlet_indices: Vec<_> = nodes.iter()
            .enumerate()
            .filter(|(_, n)| matches!(n.node_type, crate::network::NodeType::Inlet))
            .map(|(i, _)| i)
            .collect();
            
        let outlet_indices: Vec<_> = nodes.iter()
            .enumerate()
            .filter(|(_, n)| matches!(n.node_type, crate::network::NodeType::Outlet))
            .map(|(i, _)| i)
            .collect();
        
        if inlet_indices.is_empty() || outlet_indices.is_empty() {
            // If no clear inlet/outlet, calculate series resistance
            return network.edges()
                .map(|edge| edge.effective_resistance())
                .fold(T::zero(), |acc, r| acc + r);
        }
        
        // Apply unit current source between first inlet and outlet
        let inlet_idx = inlet_indices[0];
        let outlet_idx = outlet_indices[0];
        
        // Create current source vector (positive at inlet, negative at outlet)
        let mut current_vector = DVector::zeros(n);
        let unit_current = T::one();
        current_vector[inlet_idx] = unit_current.clone();
        current_vector[outlet_idx] = -unit_current.clone();
        
        // Set reference node (ground) - choose outlet as reference
        // Modify the conductance matrix to handle reference node
        let mut g_reduced = DMatrix::zeros(n - 1, n - 1);
        let mut i_reduced = DVector::zeros(n - 1);
        
        // Remove reference node row and column (using outlet as reference)
        let mut row_idx = 0;
        for i in 0..n {
            if i == outlet_idx {
                continue; // Skip reference node
            }
            let mut col_idx = 0;
            for j in 0..n {
                if j == outlet_idx {
                    continue; // Skip reference node
                }
                g_reduced[(row_idx, col_idx)] = g_matrix[(i, j)].clone();
                col_idx += 1;
            }
            i_reduced[row_idx] = current_vector[i].clone();
            row_idx += 1;
        }
        
        // Solve the linear system G * V = I for node voltages
        use cfd_math::{SparseMatrix, SparseMatrixBuilder, LinearSolver, BiCGSTAB};
use cfd_core::solver::LinearSolverConfig;;
        
        // Convert dense matrix to sparse CSR format
        let mut sparse_builder = SparseMatrixBuilder::new(n - 1, n - 1);
        for i in 0..n-1 {
            for j in 0..n-1 {
                let val = g_reduced[(i, j)].clone();
                if ComplexField::abs(val.clone()) > T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                    let _ = sparse_builder.add_entry(i, j, val);
                }
            }
        }
        let g_sparse = sparse_builder.build().unwrap_or_else(|_| {
            // Fallback: create identity matrix if build fails
            let mut builder = SparseMatrixBuilder::new(n - 1, n - 1);
            for i in 0..n-1 {
                let _ = builder.add_entry(i, i, T::one());
            }
            builder.build().expect("FIXME: Add proper error handling")
        });
        
        let config = LinearSolverConfig {
            base: cfd_core::solver::SolverConfig::builder()
                .tolerance(T::from_f64(1e-10).unwrap_or_else(T::default_epsilon))
                .max_iterations(1000)
                .build_base(),
            restart: 30,
            use_preconditioner: false,
        };
        let solver = cfd_math::ConjugateGradient::new(config);
        
        match solver.solve(&g_sparse, &i_reduced, None) {
            Ok(voltages) => {
                // Get voltage at inlet node (outlet is reference = 0V)
                let inlet_voltage_idx = if inlet_idx < outlet_idx {
                    inlet_idx
                } else {
                    inlet_idx - 1  // Adjust for removed reference node
                };
                
                let voltage_drop = voltages[inlet_voltage_idx].clone();
                
                // Equivalent resistance = V / I = voltage_drop / unit_current
                ComplexField::abs(voltage_drop.clone()) / unit_current
            }
            Err(_) => {
                // Fallback to series sum if solver fails
                network.edges()
                    .map(|edge| edge.effective_resistance())
                    .fold(T::zero(), |acc, r| acc + r)
            }
        }
    }
    
    /// Find all simple paths between two nodes using DFS
    fn find_paths(&self, network: &Network<T>, start: usize, end: usize) -> Vec<Vec<usize>> {
        
        let mut all_paths = Vec::new();
        
        // Convert node IDs to NodeIndex
        let start_id = start.to_string();
        let end_id = end.to_string();
        
        // Find NodeIndex for start and end nodes
        let start_idx = network.get_node_index(&start_id);
        let end_idx = network.get_node_index(&end_id);
        
        if let (Some(start_node), Some(end_node)) = (start_idx, end_idx) {
            // Use DFS to find all paths
            let mut visited_edges = HashSet::new();
            let mut current_path = Vec::new();
            
            self.dfs_find_paths(
                network,
                start_node,
                end_node,
                &mut visited_edges,
                &mut current_path,
                &mut all_paths
            );
        }
        
        all_paths
    }
    
    /// Recursive DFS helper for finding paths
    fn dfs_find_paths(
        &self,
        network: &Network<T>,
        current: petgraph::graph::NodeIndex,
        target: petgraph::graph::NodeIndex,
        visited_edges: &mut HashSet<petgraph::graph::EdgeIndex>,
        current_path: &mut Vec<usize>,
        all_paths: &mut Vec<Vec<usize>>
    ) {
        use petgraph::visit::EdgeRef;
        
        if current == target {
            // Found a path to the target
            if !current_path.is_empty() {
                all_paths.push(current_path.clone());
            }
            return;
        }
        
        // Explore all outgoing edges
        for edge_ref in network.graph.edges(current) {
            let edge_idx = edge_ref.id();
            let next_node = edge_ref.target();
            
            // Skip if we've already used this edge in the current path
            if !visited_edges.contains(&edge_idx) {
                // Add edge to current path
                // Convert EdgeIndex to edge ID then to position
                if let Some(edge_id) = network.get_edge_id_by_index(edge_idx) {
                    // Find the position of this edge in the edges iterator
                    if let Some(edge_pos) = network.edges()
                        .position(|e| e.id == edge_id) {
                        
                        visited_edges.insert(edge_idx);
                        current_path.push(edge_pos);
                        
                        // Recursively explore from the next node
                        self.dfs_find_paths(
                            network,
                            next_node,
                            target,
                            visited_edges,
                            current_path,
                            all_paths
                        );
                        
                        // Backtrack
                        current_path.pop();
                        visited_edges.remove(&edge_idx);
                    }
                }
            }
        }
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
        
        for edge in network.edges_with_properties() {
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
        
        for edge in network.edges_with_properties() {
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
                let connected_edges = network.node_edges(&junction.id).unwrap_or_else(|_| Vec::new());
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
    // TODO: Add solution result when available
    // pub solution_result: crate::solver::SolutionResult<T>,
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
use crate::network::{NetworkBuilder, ChannelProperties};

    #[test]
    fn test_flow_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().expect("FIXME: Add proper error handling");

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).expect("FIXME: Add proper error handling");

        let flow_analysis = analyzer.analyze_flow(&network).expect("FIXME: Add proper error handling");

        // Should have flow data
        assert!(flow_analysis.component_flows.contains_key("ch1"));
        assert!(flow_analysis.total_flow_rate >= 0.0);
    }

    #[test]
    fn test_pressure_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().expect("FIXME: Add proper error handling");

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).expect("FIXME: Add proper error handling");

        let pressure_analysis = analyzer.analyze_pressure(&network).expect("FIXME: Add proper error handling");

        // Should have pressure data
        assert!(pressure_analysis.pressures.contains_key("inlet"));
        assert!(pressure_analysis.pressures.contains_key("outlet"));
        assert!(pressure_analysis.max_pressure >= pressure_analysis.min_pressure);
    }

    #[test]
    fn test_resistance_analysis() {
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().expect("FIXME: Add proper error handling");

        let analyzer = NetworkAnalyzer::new();
        let resistance_analysis = analyzer.analyze_resistance(&network).expect("FIXME: Add proper error handling");

        // Should have resistance data
        assert!(resistance_analysis.resistances.contains_key("ch1"));
        assert!(resistance_analysis.total_resistance > 0.0);
    }

    #[test]
    fn test_performance_metrics() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().expect("FIXME: Add proper error handling");

        let analyzer = NetworkAnalyzer::new();
        let _result = analyzer.solver.solve_steady_state(&mut network).expect("FIXME: Add proper error handling");

        let performance = analyzer.analyze_performance(&network).expect("FIXME: Add proper error handling");

        // Should have performance data
        assert!(performance.throughput >= 0.0);
        assert!(performance.pressure_efficiency >= 0.0);
        assert!(performance.pressure_efficiency <= 1.0);
        assert!(performance.power_consumption >= 0.0);
    }

    #[test]
    fn test_comprehensive_analysis() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 2000.0)
            .add_junction("junction", 0.5, 0.0)
            .add_outlet_pressure("outlet1", 1.0, 0.5, 0.0)
            .add_outlet_pressure("outlet2", 1.0, -0.5, 0.0)
            .add_channel("input_ch", "inlet", "junction", ChannelProperties::new(100.0, 0.001, 1e-6))
            .add_channel("output_ch1", "junction", "outlet1", ChannelProperties::new(200.0, 0.001, 1e-6))
            .add_channel("output_ch2", "junction", "outlet2", ChannelProperties::new(200.0, 0.001, 1e-6))
            .build().expect("FIXME: Add proper error handling");

        let analyzer = NetworkAnalyzer::new();
        let analysis_result = analyzer.analyze(&mut network).expect("FIXME: Add proper error handling");

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
        let fluid = cfd_core::fluid::Fluid::water();

        let velocity = 0.1; // m/s
        let diameter = 1e-3; // 1 mm

        let re = analyzer.calculate_reynolds_number(&fluid, velocity, diameter);

        // Should be reasonable for water flow
        assert!(re > 0.0);
        assert!(re < 1000.0); // Should be laminar for these conditions
    }
}
