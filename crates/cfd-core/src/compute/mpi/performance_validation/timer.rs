//! Timing utilities for performance measurement.

use std::time::{Duration, Instant};

/// Timing utilities for performance measurement
pub struct PerformanceTimer {
    start_time: Option<Instant>,
    total_time: Duration,
}

impl PerformanceTimer {
    /// Create new performance timer
    pub fn new() -> Self {
        Self {
            start_time: None,
            total_time: Duration::default(),
        }
    }

    /// Start timing
    pub fn start(&mut self) {
        self.start_time = Some(Instant::now());
    }

    /// Stop timing and accumulate
    pub fn stop(&mut self) {
        if let Some(start) = self.start_time.take() {
            self.total_time += start.elapsed();
        }
    }

    /// Get total accumulated time
    pub fn total_time(&self) -> Duration {
        self.total_time
    }

    /// Reset timer
    pub fn reset(&mut self) {
        self.start_time = None;
        self.total_time = Duration::default();
    }
}

impl Default for PerformanceTimer {
    fn default() -> Self {
        Self::new()
    }
}
