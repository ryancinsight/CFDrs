//! Diagnostics and optimization for parallel MPI operations.
//!
//! This module provides load balancing and communication optimization
//! components for distributed CFD simulations.

#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::MpiCommunicator;
use nalgebra::RealField;

use super::LoadBalancingStrategy;

/// Load balancer for parallel matrix-free operations.
///
/// This component monitors computational load across MPI processes
/// and implements dynamic load balancing strategies.
#[cfg(feature = "mpi")]
pub struct ParallelLoadBalancer<T: RealField + Copy> {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Load balancing strategy
    strategy: LoadBalancingStrategy,
    /// Work estimates per process
    work_estimates: Vec<f64>,
    /// Communication overhead estimates
    comm_overhead: Vec<f64>,
}

#[cfg(feature = "mpi")]
impl<T: RealField + Copy> ParallelLoadBalancer<T> {
    /// Create new load balancer.
    pub fn new(communicator: &MpiCommunicator, strategy: LoadBalancingStrategy) -> Self {
        let num_procs = communicator.size() as usize;
        Self {
            communicator: communicator.clone(),
            strategy,
            work_estimates: vec![1.0; num_procs], // Start with equal load
            comm_overhead: vec![0.0; num_procs],
        }
    }

    /// Update load estimates based on recent performance.
    pub fn update_load_estimates(&mut self, local_work: f64, local_comm: f64) {
        let rank = self.communicator.rank() as usize;
        self.work_estimates[rank] = local_work;
        self.comm_overhead[rank] = local_comm;

        // Static partitioning — global communication of load/comm estimates
        // for adaptive rebalancing is deferred to the MPI-3 persistent-request path.
    }

    /// Get load balancing recommendations.
    pub fn get_recommendations(&self) -> LoadBalancingRecommendations {
        match self.strategy {
            LoadBalancingStrategy::Static => LoadBalancingRecommendations::Static,
            LoadBalancingStrategy::Dynamic => {
                // Simple dynamic balancing based on work estimates
                let avg_work: f64 =
                    self.work_estimates.iter().sum::<f64>() / self.work_estimates.len() as f64;
                let max_work = self.work_estimates.iter().cloned().fold(0.0, f64::max);
                let imbalance_ratio = max_work / avg_work;

                if imbalance_ratio > 1.2 {
                    LoadBalancingRecommendations::Rebalance
                } else {
                    LoadBalancingRecommendations::Maintain
                }
            }
            LoadBalancingStrategy::Adaptive => {
                // More sophisticated adaptive strategy
                let total_work: f64 = self.work_estimates.iter().sum();
                let total_comm: f64 = self.comm_overhead.iter().sum();
                let efficiency = total_work / (total_work + total_comm);

                if efficiency < 0.7 {
                    LoadBalancingRecommendations::OptimizeCommunication
                } else {
                    LoadBalancingRecommendations::Maintain
                }
            }
        }
    }
}

/// Load balancing recommendations.
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy)]
pub enum LoadBalancingRecommendations {
    /// Maintain current distribution
    Maintain,
    /// Rebalance load across processes
    Rebalance,
    /// Optimize communication patterns
    OptimizeCommunication,
    /// Use static distribution
    Static,
}

/// Communication optimizer for MPI operations.
///
/// This component analyzes communication patterns and suggests
/// optimizations for better parallel scaling.
#[cfg(feature = "mpi")]
pub struct CommunicationOptimizer {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Communication statistics
    stats: CommunicationStats,
}

#[cfg(feature = "mpi")]
impl CommunicationOptimizer {
    /// Create new communication optimizer.
    pub fn new(communicator: &MpiCommunicator) -> Self {
        Self {
            communicator: communicator.clone(),
            stats: CommunicationStats::default(),
        }
    }

    /// Analyze communication pattern and suggest optimizations.
    pub fn analyze_and_optimize(
        &mut self,
        message_sizes: &[usize],
        frequencies: &[f64],
    ) -> Vec<CommunicationOptimization> {
        let mut optimizations = Vec::new();

        // Analyze message size distribution
        let avg_size = message_sizes.iter().sum::<usize>() as f64 / message_sizes.len() as f64;
        let large_messages = message_sizes
            .iter()
            .filter(|&&size| size > avg_size as usize * 2)
            .count();

        if large_messages > message_sizes.len() / 4 {
            optimizations.push(CommunicationOptimization::UseNonBlockingCommunication);
        }

        // Analyze frequency patterns
        let total_freq: f64 = frequencies.iter().sum();
        let avg_freq = total_freq / frequencies.len() as f64;
        let bursty_comm = frequencies
            .iter()
            .filter(|&&freq| freq > avg_freq * 2.0)
            .count();

        if bursty_comm > frequencies.len() / 3 {
            optimizations.push(CommunicationOptimization::ImplementCommunicationAggregation);
        }

        // Check for potential overlap opportunities
        if self.stats.computation_time > self.stats.communication_time * 2.0 {
            optimizations.push(CommunicationOptimization::EnableComputationCommunicationOverlap);
        }

        optimizations
    }

    /// Update communication statistics.
    pub fn update_stats(&mut self, comp_time: f64, comm_time: f64, messages_sent: usize) {
        self.stats.computation_time = comp_time;
        self.stats.communication_time = comm_time;
        self.stats.messages_sent = messages_sent;
    }
}

/// Communication statistics.
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Default)]
struct CommunicationStats {
    /// Time spent in computation
    computation_time: f64,
    /// Time spent in communication
    communication_time: f64,
    /// Number of messages sent
    messages_sent: usize,
}

/// Communication optimization recommendations.
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy)]
pub enum CommunicationOptimization {
    /// Use non-blocking MPI operations
    UseNonBlockingCommunication,
    /// Aggregate small messages
    ImplementCommunicationAggregation,
    /// Overlap computation with communication
    EnableComputationCommunicationOverlap,
    /// Use one-sided communication
    UseOneSidedCommunication,
}
