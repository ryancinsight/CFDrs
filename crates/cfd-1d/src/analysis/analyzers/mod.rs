//! Domain-specific network analyzers

mod flow;
mod performance;
mod pressure;
mod resistance;
mod traits;

pub use flow::FlowAnalyzer;
pub use performance::PerformanceAnalyzer;
pub use pressure::PressureAnalyzer;
pub use resistance::ResistanceAnalyzer;
pub use traits::NetworkAnalyzer;

use super::{FlowAnalysis, PerformanceMetrics, PressureAnalysis, ResistanceAnalysis};
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
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

/// Comprehensive network analyzer orchestrator
pub struct NetworkAnalyzerOrchestrator<T: RealField + Copy> {
    flow_analyzer: FlowAnalyzer<T>,
    pressure_analyzer: PressureAnalyzer<T>,
    resistance_analyzer: ResistanceAnalyzer<T>,
    performance_analyzer: PerformanceAnalyzer<T>,
    solver: crate::solver::NetworkSolver<T>,
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> Default for NetworkAnalyzerOrchestrator<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzerOrchestrator<T> {
    /// Create a new network analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            flow_analyzer: FlowAnalyzer::new(),
            pressure_analyzer: PressureAnalyzer::new(),
            resistance_analyzer: ResistanceAnalyzer::new(),
            performance_analyzer: PerformanceAnalyzer::new(),
            solver: crate::solver::NetworkSolver::new(),
        }
    }

    /// Create analyzer with custom solver configuration
    pub fn with_solver_config(config: crate::solver::SolverConfig<T>) -> Self {
        Self {
            flow_analyzer: FlowAnalyzer::new(),
            pressure_analyzer: PressureAnalyzer::new(),
            resistance_analyzer: ResistanceAnalyzer::new(),
            performance_analyzer: PerformanceAnalyzer::new(),
            solver: crate::solver::NetworkSolver::with_config(config),
        }
    }

    /// Perform comprehensive network analysis
    pub fn analyze(
        &mut self,
        network: &mut crate::network::Network<T>,
    ) -> Result<NetworkAnalysisResult<T>> {
        // Solve the network
        let problem = crate::solver::NetworkProblem::new(network.clone());
        let solved_network = self.solver.solve_network(&problem)?;

        // Perform domain-specific analyses
        let flow_analysis = self.flow_analyzer.analyze(&solved_network)?;
        let pressure_analysis = self.pressure_analyzer.analyze(&solved_network)?;
        let resistance_analysis = self.resistance_analyzer.analyze(&solved_network)?;
        let performance_metrics = self.performance_analyzer.analyze(&solved_network)?;

        // Update the original network with solution
        *network = solved_network;

        Ok(NetworkAnalysisResult {
            flow_analysis,
            pressure_analysis,
            resistance_analysis,
            performance_metrics,
        })
    }
}
