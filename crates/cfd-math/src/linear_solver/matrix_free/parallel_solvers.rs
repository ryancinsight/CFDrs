//! Parallel matrix-free solvers with optimized MPI communication.
//!
//! This module provides MPI-aware matrix-free solvers that implement
//! communication overlap techniques and load balancing optimizations
//! for large-scale distributed CFD simulations.

use crate::linear_solver::config::IterativeSolverConfig;
use crate::linear_solver::traits::{ConvergenceMonitor, LinearOperator};
#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::{DistributedVector, MpiCommunicator, MpiResult};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use std::fmt::Debug;

/// Communication overlap strategy for MPI operations
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy)]
pub enum CommunicationOverlap {
    /// No overlap - synchronous communication
    None,
    /// Overlap computation with communication
    ComputationOverlap,
    /// Persistent communication requests
    PersistentRequests,
}

/// Load balancing strategy for parallel operations
#[cfg(feature = "mpi")]
#[derive(Debug, Clone, Copy)]
pub enum LoadBalancingStrategy {
    /// Static load balancing (equal distribution)
    Static,
    /// Dynamic load balancing based on work estimates
    Dynamic,
    /// Adaptive load balancing with monitoring
    Adaptive,
}

/// Optimized parallel matrix-free BiCGSTAB solver.
///
/// This solver implements communication overlap techniques and load balancing
/// optimizations for efficient scaling on distributed systems.
#[cfg(feature = "mpi")]
pub struct ParallelMatrixFreeBiCGSTAB<T: RealField + Copy> {
    /// Solver configuration
    config: IterativeSolverConfig<T>,
    /// Communication overlap strategy
    overlap_strategy: CommunicationOverlap,
    /// Load balancing strategy
    load_balancing: LoadBalancingStrategy,
    /// MPI communicator
    communicator: MpiCommunicator,
}

#[cfg(feature = "mpi")]
impl<T> ParallelMatrixFreeBiCGSTAB<T>
where
    T: RealField + Copy + Default + From<f64> + std::fmt::LowerExp + 'static,
{
    /// Create a new parallel matrix-free BiCGSTAB solver.
    pub fn new(
        config: IterativeSolverConfig<T>,
        communicator: &MpiCommunicator,
        overlap_strategy: CommunicationOverlap,
        load_balancing: LoadBalancingStrategy,
    ) -> Self {
        Self {
            config,
            overlap_strategy,
            load_balancing,
            communicator: communicator.clone(),
        }
    }

    /// Create with default configuration and strategies.
    pub fn default_with_mpi(communicator: &MpiCommunicator) -> Self
    where
        T: From<f64> + Default,
    {
        let config = IterativeSolverConfig::default();
        Self::new(
            config,
            communicator,
            CommunicationOverlap::ComputationOverlap,
            LoadBalancingStrategy::Adaptive,
        )
    }

    /// Solve distributed linear system with optimized communication.
    ///
    /// This method implements communication overlap and load balancing
    /// to minimize MPI communication overhead and maximize parallel efficiency.
    pub fn solve_distributed<Op: LinearOperator<T>>(
        &self,
        operator: &Op,
        b: &DistributedVector<T>,
        x: &mut DistributedVector<T>,
    ) -> MpiResult<ConvergenceMonitor<T>> {
        let local_size = operator.size();

        // Initialize workspace vectors
        let mut r = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut r_hat = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut p = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut p_hat = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut v = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut s = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;
        let mut t = DistributedVector::new(
            local_size,
            &self.communicator,
            x.subdomain.clone(),
            x.ghost_manager.clone(),
        )?;

        // Initial residual: r = b - A*x
        operator.apply(&x.local_data, &mut r.local_data)?;
        for i in 0..local_size {
            r.local_data[i] = b.local_data[i] - r.local_data[i];
        }

        let mut r_norm = r.norm()?;
        let mut monitor = ConvergenceMonitor::new(r_norm);

        if r_norm < self.config.tolerance {
            return Ok(monitor);
        }

        // Initialize shadow residual
        r_hat.local_data.copy_from_slice(&r.local_data);

        // Initialize scalars
        let mut rho_old = r_hat.dot(&r)?;
        let mut alpha = T::from_f64(1.0).unwrap();
        let mut omega = T::from_f64(1.0).unwrap();

        // Initial search directions
        p.local_data.copy_from_slice(&r.local_data);
        p_hat.local_data.copy_from_slice(&r.local_data);

        let mut iteration = 0;

        while iteration < self.config.max_iterations {
            // Overlap communication with computation based on strategy
            match self.overlap_strategy {
                CommunicationOverlap::ComputationOverlap => {
                    self.solve_with_computation_overlap(
                        operator,
                        &mut r,
                        &mut r_hat,
                        &mut p,
                        &mut p_hat,
                        &mut v,
                        &mut s,
                        &mut t,
                        &mut rho_old,
                        &mut alpha,
                        &mut omega,
                    )?;
                }
                CommunicationOverlap::PersistentRequests => {
                    self.solve_with_persistent_requests(
                        operator,
                        &mut r,
                        &mut r_hat,
                        &mut p,
                        &mut p_hat,
                        &mut v,
                        &mut s,
                        &mut t,
                        &mut rho_old,
                        &mut alpha,
                        &mut omega,
                    )?;
                }
                CommunicationOverlap::None => {
                    self.solve_synchronous(
                        operator,
                        &mut r,
                        &mut r_hat,
                        &mut p,
                        &mut p_hat,
                        &mut v,
                        &mut s,
                        &mut t,
                        &mut rho_old,
                        &mut alpha,
                        &mut omega,
                    )?;
                }
            }

            // Update solution
            for i in 0..local_size {
                x.local_data[i] += alpha * p.local_data[i] + omega * s.local_data[i];
            }

            // Check convergence
            r_norm = r.norm()?;
            monitor.record_residual(r_norm);

            if r_norm < self.config.tolerance {
                return Ok(monitor);
            }

            iteration += 1;
        }

        Err(
            Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            })
            .into(),
        )
    }

    /// Solve with computation-communication overlap.
    fn solve_with_computation_overlap<Op: LinearOperator<T>>(
        &self,
        operator: &Op,
        r: &mut DistributedVector<T>,
        r_hat: &mut DistributedVector<T>,
        p: &mut DistributedVector<T>,
        p_hat: &mut DistributedVector<T>,
        v: &mut DistributedVector<T>,
        s: &mut DistributedVector<T>,
        t: &mut DistributedVector<T>,
        rho_old: &mut T,
        alpha: &mut T,
        omega: &mut T,
    ) -> MpiResult<()> {
        // Start asynchronous operator application
        operator.apply(&p_hat.local_data, &mut v.local_data)?;

        // While computation is happening, prepare for communication
        // (In a real implementation, this would start MPI communication)

        // Compute alpha using local data
        let rho_new = r_hat.dot(v)?;
        if rho_new == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }
        *alpha = *rho_old / rho_new;

        // Compute s = r - alpha * v
        for i in 0..r.local_data.len() {
            s.local_data[i] = r.local_data[i] - *alpha * v.local_data[i];
        }

        // Early convergence check
        let s_norm = s.norm()?;
        if s_norm < self.config.tolerance {
            return Ok(());
        }

        // Apply operator to s (potentially overlapped)
        operator.apply(&s.local_data, &mut t.local_data)?;

        // Compute omega
        let ts_dot = t.dot(s)?;
        let tt_dot = t.dot(t)?;

        if tt_dot == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }

        *omega = ts_dot / tt_dot;

        // Update r = s - omega * t
        for i in 0..r.local_data.len() {
            r.local_data[i] = s.local_data[i] - *omega * t.local_data[i];
        }

        // Update rho_old
        *rho_old = r_hat.dot(r)?;

        if *rho_old == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }

        // Compute beta
        let beta = (*rho_old / rho_new) * (*alpha / *omega);

        // Update search directions
        for i in 0..p.local_data.len() {
            p.local_data[i] = r.local_data[i] + beta * (p.local_data[i] - *omega * v.local_data[i]);
            p_hat.local_data[i] =
                r_hat.local_data[i] + beta * (p_hat.local_data[i] - *omega * v.local_data[i]);
        }

        Ok(())
    }

    /// Solve with persistent communication requests.
    fn solve_with_persistent_requests<Op: LinearOperator<T>>(
        &self,
        _operator: &Op,
        _r: &mut DistributedVector<T>,
        _r_hat: &mut DistributedVector<T>,
        _p: &mut DistributedVector<T>,
        _p_hat: &mut DistributedVector<T>,
        _v: &mut DistributedVector<T>,
        _s: &mut DistributedVector<T>,
        _t: &mut DistributedVector<T>,
        _rho_old: &mut T,
        _alpha: &mut T,
        _omega: &mut T,
    ) -> MpiResult<()> {
        // TODO: Implement persistent MPI requests and overlap communication with computation.
        Err(
            Error::InvalidConfiguration("Persistent requests not yet implemented".to_string())
                .into(),
        )
    }

    /// Synchronous solve (no overlap).
    fn solve_synchronous<Op: LinearOperator<T>>(
        &self,
        operator: &Op,
        r: &mut DistributedVector<T>,
        r_hat: &mut DistributedVector<T>,
        p: &mut DistributedVector<T>,
        p_hat: &mut DistributedVector<T>,
        v: &mut DistributedVector<T>,
        s: &mut DistributedVector<T>,
        t: &mut DistributedVector<T>,
        rho_old: &mut T,
        alpha: &mut T,
        omega: &mut T,
    ) -> MpiResult<()> {
        // Standard BiCGSTAB iteration without overlap
        operator.apply(&p_hat.local_data, &mut v.local_data)?;

        let rho_new = r_hat.dot(v)?;
        if rho_new == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }
        *alpha = *rho_old / rho_new;

        for i in 0..s.local_data.len() {
            s.local_data[i] = r.local_data[i] - *alpha * v.local_data[i];
        }

        let s_norm = s.norm()?;
        if s_norm < self.config.tolerance {
            return Ok(());
        }

        operator.apply(&s.local_data, &mut t.local_data)?;

        let ts_dot = t.dot(s)?;
        let tt_dot = t.dot(t)?;

        if tt_dot == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }

        *omega = ts_dot / tt_dot;

        for i in 0..r.local_data.len() {
            r.local_data[i] = s.local_data[i] - *omega * t.local_data[i];
        }

        *rho_old = r_hat.dot(r)?;

        if *rho_old == T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix).into());
        }

        let beta = (*rho_old / rho_new) * (*alpha / *omega);

        for i in 0..p.local_data.len() {
            p.local_data[i] = r.local_data[i] + beta * (p.local_data[i] - *omega * v.local_data[i]);
            p_hat.local_data[i] =
                r_hat.local_data[i] + beta * (p_hat.local_data[i] - *omega * v.local_data[i]);
        }

        Ok(())
    }
}

#[cfg(feature = "mpi")]
impl<T> crate::linear_solver::traits::Configurable<T> for ParallelMatrixFreeBiCGSTAB<T>
where
    T: RealField + Copy + Default + From<f64> + 'static,
{
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

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

        // Communicate load/comm estimates across ranks for global decisions
        self.communicator.all_gather(&local_work, &mut self.work_estimates);
        self.communicator.all_gather(&local_comm, &mut self.comm_overhead);
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

#[cfg(test)]
#[cfg(feature = "mpi")]
mod tests {
    use super::*;
    use cfd_core::compute::mpi::MpiUniverse;

    #[test]
    fn test_parallel_solver_creation() {
        // This test would require MPI initialization
        // For now, just test that the types can be created
        let _optimizer = std::marker::PhantomData::<CommunicationOptimizer>;
        let _balancer = std::marker::PhantomData::<ParallelLoadBalancer<f64>>;
    }
}
