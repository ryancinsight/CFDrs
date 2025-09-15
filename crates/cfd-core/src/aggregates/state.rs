//! Simulation state management

use serde::{Deserialize, Serialize};

/// Simulation state enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimulationState {
    /// Initial state
    Initialized,
    /// Simulation is running
    Running,
    /// Simulation is paused
    Paused,
    /// Simulation completed successfully
    Completed,
    /// Simulation failed with error
    Failed,
}

impl Default for SimulationState {
    fn default() -> Self {
        Self::Initialized
    }
}

impl SimulationState {
    /// Check if simulation is active
    pub fn is_active(&self) -> bool {
        matches!(self, Self::Running | Self::Paused)
    }

    /// Check if simulation has finished
    pub fn is_finished(&self) -> bool {
        matches!(self, Self::Completed | Self::Failed)
    }

    /// Transition to running state
    ///
    /// # Errors
    /// Returns an error if:
    /// - Current state is not Initialized or Paused
    pub fn start(&mut self) -> Result<(), &'static str> {
        match self {
            Self::Initialized | Self::Paused => {
                *self = Self::Running;
                Ok(())
            }
            _ => Err("Cannot start simulation from current state"),
        }
    }

    /// Transition to paused state
    pub fn pause(&mut self) -> Result<(), &'static str> {
        match self {
            Self::Running => {
                *self = Self::Paused;
                Ok(())
            }
            _ => Err("Can only pause running simulation"),
        }
    }

    /// Transition to completed state
    pub fn complete(&mut self) -> Result<(), &'static str> {
        match self {
            Self::Running => {
                *self = Self::Completed;
                Ok(())
            }
            _ => Err("Can only complete running simulation"),
        }
    }

    /// Transition to failed state
    pub fn fail(&mut self) {
        *self = Self::Failed;
    }
}
