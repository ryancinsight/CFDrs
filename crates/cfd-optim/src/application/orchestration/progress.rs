//! Thread-safe progress tracker for parallel candidate evaluation.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Logs periodic heartbeat messages during parallel evaluation scans.
///
/// Use with `rayon::par_iter()`:
/// ```rust,ignore
/// let progress = Arc::new(ScanProgress::new("option2 scan", candidates.len()));
/// let progress_ref = Arc::clone(&progress);
/// candidates.par_iter().for_each(move |c| {
///     // ... evaluate ...
///     progress_ref.record();
/// });
/// progress.finish();
/// ```
pub struct ScanProgress {
    label: &'static str,
    total: usize,
    started_at: Instant,
    processed: AtomicUsize,
    next_log_at: AtomicUsize,
    log_step: usize,
}

impl ScanProgress {
    /// Create and log the initial "starting N evaluations" message.
    pub fn new(label: &'static str, total: usize) -> Self {
        let log_step = if total <= 25 {
            1
        } else {
            (total / 20).clamp(25, 1_000)
        };
        tracing::info!(
            "      {label}: starting {} evaluations (heartbeat every {})",
            total,
            log_step
        );
        Self {
            label,
            total,
            started_at: Instant::now(),
            processed: AtomicUsize::new(0),
            next_log_at: AtomicUsize::new(log_step.min(total.max(1))),
            log_step,
        }
    }

    /// Record one completed evaluation; log if the heartbeat threshold is reached.
    pub fn record(&self) {
        let processed = self.processed.fetch_add(1, Ordering::Relaxed) + 1;
        loop {
            let next_log_at = self.next_log_at.load(Ordering::Relaxed);
            if processed < next_log_at && processed < self.total {
                break;
            }
            let new_next = (next_log_at + self.log_step).min(self.total.max(1));
            if self
                .next_log_at
                .compare_exchange(next_log_at, new_next, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                self.log(processed);
                break;
            }
        }
    }

    /// Log the final count after the scan completes.
    pub fn finish(&self) {
        self.log(self.processed.load(Ordering::Relaxed));
    }

    fn log(&self, processed: usize) {
        let elapsed = self.started_at.elapsed().as_secs_f64().max(1.0e-9);
        let pct = if self.total == 0 {
            100.0
        } else {
            100.0 * processed as f64 / self.total as f64
        };
        let rate = processed as f64 / elapsed;
        tracing::info!(
            "      {}: {}/{} ({:.1}%) in {:.1}s [{:.1} candidates/s]",
            self.label,
            processed,
            self.total,
            pct,
            elapsed,
            rate
        );
    }
}
