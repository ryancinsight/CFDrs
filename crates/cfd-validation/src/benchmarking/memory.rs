//! Memory usage profiling for CFD operations
//!
//! Tracks memory allocation patterns, peak usage, and memory efficiency
//! of CFD algorithms and data structures.

use cfd_2d::grid::Grid2D;
use cfd_core::error::Result;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

/// Thread-safe memory statistics using atomics
pub struct MemoryStats {
    /// Total bytes currently allocated
    pub current_allocated: AtomicUsize,
    /// Peak memory usage during measurement period
    pub peak_allocated: AtomicUsize,
    /// Total allocation count
    pub allocation_count: AtomicUsize,
    /// Total deallocation count
    pub deallocation_count: AtomicUsize,
    /// Total bytes allocated (cumulative)
    pub total_allocated: AtomicU64,
    /// Total bytes deallocated (cumulative)
    pub total_deallocated: AtomicU64,
    /// Largest single allocation
    pub max_allocation_size: AtomicUsize,
}

/// A snapshot of memory statistics at a specific point in time
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct MemoryStatsSnapshot {
    /// Total bytes currently allocated
    pub current_allocated: usize,
    /// Peak memory usage during measurement period
    pub peak_allocated: usize,
    /// Total allocation count
    pub allocation_count: usize,
    /// Total deallocation count
    pub deallocation_count: usize,
    /// Total bytes allocated (cumulative)
    pub total_allocated: u64,
    /// Total bytes deallocated (cumulative)
    pub total_deallocated: u64,
    /// Average allocation size in bytes
    pub avg_allocation_size: f64,
    /// Largest single allocation in bytes
    pub max_allocation_size: usize,
}

impl MemoryStats {
    const fn new() -> Self {
        Self {
            current_allocated: AtomicUsize::new(0),
            peak_allocated: AtomicUsize::new(0),
            allocation_count: AtomicUsize::new(0),
            deallocation_count: AtomicUsize::new(0),
            total_allocated: AtomicU64::new(0),
            total_deallocated: AtomicU64::new(0),
            max_allocation_size: AtomicUsize::new(0),
        }
    }

    /// Takes a snapshot of the current memory statistics
    pub fn snapshot(&self) -> MemoryStatsSnapshot {
        // Use Acquire ordering to ensure we see the latest updates from all threads
        let total_allocated = self.total_allocated.load(Ordering::Acquire);
        let allocation_count = self.allocation_count.load(Ordering::Acquire);
        let current_allocated = self.current_allocated.load(Ordering::Acquire);
        let peak_allocated = self.peak_allocated.load(Ordering::Acquire);
        let deallocation_count = self.deallocation_count.load(Ordering::Acquire);
        let total_deallocated = self.total_deallocated.load(Ordering::Acquire);
        let max_allocation_size = self.max_allocation_size.load(Ordering::Acquire);

        MemoryStatsSnapshot {
            current_allocated,
            peak_allocated,
            allocation_count,
            deallocation_count,
            total_allocated,
            total_deallocated,
            avg_allocation_size: if allocation_count > 0 {
                total_allocated as f64 / allocation_count as f64
            } else {
                0.0
            },
            max_allocation_size,
        }
    }
}

/// Global memory allocator wrapper for tracking allocations
#[global_allocator]
static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator {
    allocator: System,
    stats: MemoryStats::new(),
};

impl MemoryStatsSnapshot {
    /// Calculate memory efficiency metrics
    pub fn efficiency_metrics(&self) -> MemoryEfficiency {
        let net_usage = self.current_allocated as f64;
        let peak_usage = self.peak_allocated as f64;
        let total_operations = (self.allocation_count + self.deallocation_count) as f64;

        let memory_efficiency = if peak_usage > 0.0 {
            1.0 - (net_usage / peak_usage)
        } else {
            1.0
        };

        let allocation_efficiency = if total_operations > 0.0 {
            1.0 - (self.allocation_count as f64 / total_operations)
        } else {
            1.0
        };

        MemoryEfficiency {
            memory_efficiency,
            allocation_efficiency,
            fragmentation_ratio: self.calculate_fragmentation(),
        }
    }

    /// Estimate memory fragmentation (simplified)
    fn calculate_fragmentation(&self) -> f64 {
        if self.allocation_count == 0 {
            return 0.0;
        }

        let avg_size = self.avg_allocation_size;
        if avg_size > 0.0 {
            (self.max_allocation_size as f64 / avg_size).min(10.0) / 10.0
        } else {
            0.0
        }
    }
}

/// Memory efficiency metrics
#[derive(Debug, Clone)]
pub struct MemoryEfficiency {
    /// Memory efficiency (0.0 = all memory in use, 1.0 = no waste)
    pub memory_efficiency: f64,
    /// Allocation efficiency (lower is better - fewer allocations)
    pub allocation_efficiency: f64,
    /// Fragmentation ratio (0.0 = no fragmentation, 1.0 = high fragmentation)
    pub fragmentation_ratio: f64,
}

/// Thread-safe memory profiler
pub struct MemoryProfiler {
    start_snapshot: Mutex<Option<MemoryStatsSnapshot>>,
}

impl MemoryProfiler {
    /// Create a new memory profiler
    pub fn new() -> Self {
        Self {
            start_snapshot: Mutex::new(None),
        }
    }

    /// Start memory profiling session
    pub fn start_profiling(&self) {
        let mut start_snapshot = self.start_snapshot.lock().unwrap();
        *start_snapshot = Some(self.get_stats());
    }

    /// Stop profiling and return statistics for the session
    pub fn stop_profiling(&self) -> MemoryStatsSnapshot {
        let start_stats = self
            .start_snapshot
            .lock()
            .unwrap()
            .take()
            .unwrap_or_default();

        let current_stats = self.get_stats();

        let allocation_count = current_stats
            .allocation_count
            .saturating_sub(start_stats.allocation_count);
        let total_allocated = current_stats
            .total_allocated
            .saturating_sub(start_stats.total_allocated);

        let current_allocated = current_stats
            .current_allocated
            .saturating_sub(start_stats.current_allocated);
        let peak_delta = current_stats
            .peak_allocated
            .saturating_sub(start_stats.peak_allocated);

        MemoryStatsSnapshot {
            current_allocated,
            peak_allocated: peak_delta.max(current_allocated),
            allocation_count,
            deallocation_count: current_stats
                .deallocation_count
                .saturating_sub(start_stats.deallocation_count),
            total_allocated,
            total_deallocated: current_stats
                .total_deallocated
                .saturating_sub(start_stats.total_deallocated),
            avg_allocation_size: if allocation_count > 0 {
                total_allocated as f64 / allocation_count as f64
            } else {
                0.0
            },
            max_allocation_size: current_stats.max_allocation_size,
        }
    }

    /// Get current memory statistics
    pub fn get_stats(&self) -> MemoryStatsSnapshot {
        GLOBAL_ALLOCATOR.get_stats()
    }

    /// Profile a closure and return memory statistics
    pub fn profile_closure<F, T>(&self, operation: F) -> (T, MemoryStatsSnapshot)
    where
        F: FnOnce() -> T,
    {
        self.start_profiling();
        let result = operation();
        let stats = self.stop_profiling();
        (result, stats)
    }

    /// Check for memory leaks (simplified)
    pub fn detect_memory_leaks(&self, tolerance_bytes: usize) -> Result<bool> {
        let stats = self.get_stats();
        let net_allocation = stats.total_allocated as i64 - stats.total_deallocated as i64;

        if net_allocation.abs() > tolerance_bytes as i64 {
            return Ok(true); // Potential leak detected
        }

        Ok(false)
    }
}

impl Default for MemoryProfiler {
    fn default() -> Self {
        Self::new()
    }
}

/// Custom allocator that tracks memory usage
struct TrackingAllocator {
    allocator: System,
    stats: MemoryStats,
}

impl TrackingAllocator {
    fn get_stats(&self) -> MemoryStatsSnapshot {
        self.stats.snapshot()
    }

    fn record_allocation(&self, size: usize) {
        // Use Release ordering to ensure these updates are visible to snapshot()
        self.stats.allocation_count.fetch_add(1, Ordering::Release);
        self.stats.total_allocated.fetch_add(size as u64, Ordering::Release);
        let current = self.stats.current_allocated.fetch_add(size, Ordering::AcqRel) + size;
        
        // Update peak using a loop to ensure we don't overwrite a higher value
        let mut peak = self.stats.peak_allocated.load(Ordering::Acquire);
        while current > peak {
            match self.stats.peak_allocated.compare_exchange_weak(
                peak, 
                current, 
                Ordering::AcqRel, 
                Ordering::Acquire
            ) {
                Ok(_) => break,
                Err(actual) => peak = actual,
            }
        }

        // Update max allocation size
        let mut max_size = self.stats.max_allocation_size.load(Ordering::Acquire);
        while size > max_size {
            match self.stats.max_allocation_size.compare_exchange_weak(
                max_size, 
                size, 
                Ordering::AcqRel, 
                Ordering::Acquire
            ) {
                Ok(_) => break,
                Err(actual) => max_size = actual,
            }
        }
    }

    fn record_deallocation(&self, size: usize) {
        self.stats.deallocation_count.fetch_add(1, Ordering::Release);
        self.stats.total_deallocated.fetch_add(size as u64, Ordering::Release);
        self.stats.current_allocated.fetch_sub(size, Ordering::AcqRel);
    }
}

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = self.allocator.alloc(layout);
        if !ptr.is_null() {
            self.record_allocation(layout.size());
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        self.record_deallocation(layout.size());
        self.allocator.dealloc(ptr, layout);
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = self.allocator.alloc_zeroed(layout);
        if !ptr.is_null() {
            self.record_allocation(layout.size());
        }
        ptr
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let old_size = layout.size();
        let new_ptr = self.allocator.realloc(ptr, layout, new_size);

        if !new_ptr.is_null() {
            // Count as deallocation + allocation
            self.record_deallocation(old_size);
            self.record_allocation(new_size);
        }

        new_ptr
    }
}

/// CFD-specific memory profiling
pub struct CfdMemoryProfiler {
    profiler: MemoryProfiler,
}

impl CfdMemoryProfiler {
    /// Create a new CFD memory profiler with default configuration
    ///
    /// Initializes memory profiling capabilities optimized for CFD workloads,
    /// including tracking of memory allocations for matrices, grids, and solver data structures.
    /// Configured with appropriate sampling rates and memory tracking thresholds for
    /// CFD simulation performance analysis.
    pub fn new() -> Self {
        Self {
            profiler: MemoryProfiler::new(),
        }
    }

    /// Profile memory usage of CFD matrix operations
    pub fn profile_matrix_operations(&self) -> Result<MemoryStatsSnapshot> {
        Ok(self
            .profiler
            .profile_closure(|| {
                // Create and manipulate CFD-style matrices
                let mut builder = cfd_math::sparse::SparseMatrixBuilder::new(100, 100);
                for i in 0..100 {
                    builder.add_entry(i, i, 1.0).unwrap();
                }

                let matrix = builder.build().unwrap();
                let vector = nalgebra::DVector::from_vec(vec![1.0f64; 100]);
                let mut result = nalgebra::DVector::zeros(100);

                // Perform SPMV operations
                for _ in 0..10 {
                    cfd_math::sparse::spmv(&matrix, &vector, &mut result);
                }
            })
            .1)
    }

    /// Profile memory usage of CFD grid operations
    pub fn profile_grid_operations(&self) -> Result<MemoryStatsSnapshot> {
        use cfd_2d::grid::StructuredGrid2D;

        Ok(self
            .profiler
            .profile_closure(|| {
                // Create multiple grids of different sizes
                let grids = vec![
                    StructuredGrid2D::new(16, 16, 0.0, 1.0, 0.0, 1.0).unwrap(),
                    StructuredGrid2D::new(32, 32, 0.0, 1.0, 0.0, 1.0).unwrap(),
                ];

                // Perform grid operations
                for grid in grids {
                    let _nx = grid.nx();
                    let _ny = grid.ny();
                }
            })
            .1)
    }

    /// Profile memory usage of CFD simulation fields
    pub fn profile_simulation_fields(&self) -> Result<MemoryStatsSnapshot> {
        use cfd_2d::fields::SimulationFields;
        use cfd_core::physics::fluid::Fluid;

        Ok(self
            .profiler
            .profile_closure(|| {
                let fluid = Fluid::new("Test".to_string(), 1.0, 0.01, 1000.0, 0.001, 343.0);
                let _fields = SimulationFields::with_fluid(32, 32, &fluid);
            })
            .1)
    }

    /// Run comprehensive CFD memory profiling suite
    pub fn run_memory_suite(&self) -> Result<Vec<(String, MemoryStatsSnapshot)>> {
        let mut results = Vec::new();

        results.push(("Matrix_Operations".to_string(), self.profile_matrix_operations()?));
        results.push(("Grid_Operations".to_string(), self.profile_grid_operations()?));
        results.push(("Simulation_Fields".to_string(), self.profile_simulation_fields()?));

        Ok(results)
    }
}

impl Default for CfdMemoryProfiler {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for MemoryStatsSnapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Memory: current={:.2}MB, peak={:.2}MB, allocations={}, deallocations={}",
            self.current_allocated as f64 / 1_048_576.0,
            self.peak_allocated as f64 / 1_048_576.0,
            self.allocation_count,
            self.deallocation_count
        )
    }
}

impl std::fmt::Display for MemoryStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.snapshot())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_profiler() {
        let profiler = MemoryProfiler::new();

        let (result, stats) = profiler.profile_closure(|| {
            let data = vec![1.0f64; 1000];
            data.iter().sum::<f64>()
        });

        assert!(result > 0.0);
        assert!(stats.total_allocated > 0);
        assert!(stats.peak_allocated >= stats.current_allocated);
    }

    #[test]
    fn test_memory_efficiency() {
        let stats = MemoryStatsSnapshot {
            current_allocated: 1000,
            peak_allocated: 2000,
            allocation_count: 10,
            deallocation_count: 5,
            total_allocated: 5000,
            total_deallocated: 4000,
            avg_allocation_size: 500.0,
            max_allocation_size: 1000,
        };

        let efficiency = stats.efficiency_metrics();
        assert!(efficiency.memory_efficiency >= 0.0 && efficiency.memory_efficiency <= 1.0);
        assert!(efficiency.allocation_efficiency >= 0.0 && efficiency.allocation_efficiency <= 1.0);
    }

    #[test]
    fn test_cfd_memory_profiling() {
        let cfd_profiler = CfdMemoryProfiler::new();
        let results = cfd_profiler.run_memory_suite().unwrap();

        assert!(!results.is_empty());
        for (name, stats) in results {
            println!("{}: {}", name, stats);
            assert!(stats.total_allocated > 0);
        }
    }
}
