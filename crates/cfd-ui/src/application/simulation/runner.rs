//! Simulation runner — executes CFD simulation on a background thread.

use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::Arc;

/// State of a simulation run.
#[derive(Clone, Debug)]
pub enum SimulationState {
    /// No simulation running.
    Idle,
    /// Simulation is running.
    Running,
    /// Simulation completed successfully.
    Completed,
    /// Simulation failed with an error message.
    Failed(String),
}

/// Tracks progress and cancellation of a background simulation.
pub struct SimulationRunner {
    /// Current state.
    pub state: SimulationState,
    /// Progress percentage (0..100).
    pub progress: Arc<AtomicU32>,
    /// Cancellation flag.
    pub cancel: Arc<AtomicBool>,
}

impl SimulationRunner {
    /// Create a new idle runner.
    #[must_use]
    pub fn new() -> Self {
        Self {
            state: SimulationState::Idle,
            progress: Arc::new(AtomicU32::new(0)),
            cancel: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Get the current progress percentage.
    #[must_use]
    pub fn progress_percent(&self) -> u32 {
        self.progress.load(Ordering::Relaxed)
    }

    /// Request cancellation of the running simulation.
    pub fn request_cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// Check if cancellation was requested.
    #[must_use]
    pub fn is_cancel_requested(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// Reset the runner to idle state.
    pub fn reset(&mut self) {
        self.state = SimulationState::Idle;
        self.progress.store(0, Ordering::Relaxed);
        self.cancel.store(false, Ordering::Relaxed);
    }
}

impl Default for SimulationRunner {
    fn default() -> Self {
        Self::new()
    }
}
