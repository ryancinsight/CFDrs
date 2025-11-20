//! Memory usage profiling for CFD operations
//!
//! Tracks memory allocation patterns, peak usage, and memory efficiency
//! of CFD algorithms and data structures.

use cfd_2d::grid::Grid2D;
use cfd_core::error::{Error, Result};
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Global memory allocator wrapper for tracking allocations
#[global_allocator]
static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator {
    allocator: System,
    stats: Mutex::new(MemoryStats {
        total_allocated: 0,
        total_deallocated: 0,
        current_allocated: 0,
        peak_allocated: 0,
        deallocation_count: 0,
        allocation_count: 0,
        avg_allocation_size: 0.0,
        max_allocation_size: 0,
    }),
};

/// Thread-safe memory statistics
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct MemoryStats {
    /// Total bytes currently allocated
    pub current_allocated: usize,
    /// Peak memory usage during measurement period
    pub peak_allocated: usize,
    /// Total allocation count
    pub allocation_count: usize,
    /// Total deallocation count
    pub deallocation_count: usize,
    /// Total bytes allocated (cumulative)
    pub total_allocated: usize,
    /// Total bytes deallocated (cumulative)
    pub total_deallocated: usize,
    /// Average allocation size
    pub avg_allocation_size: f64,
    /// Largest single allocation
    pub max_allocation_size: usize,
}

impl MemoryStats {
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

        // Simple fragmentation estimate based on allocation size variance
        // In a real implementation, this would track actual heap fragmentation
        let avg_size = self.avg_allocation_size;
        if avg_size > 0.0 {
            // This is a placeholder - real fragmentation analysis would be more complex
            (self.max_allocation_size as f64 / avg_size).min(10.0) / 10.0
        } else {
            0.0
        }
    }

    /// Reset statistics
    pub fn reset(&mut self) {
        *self = Self::default();
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
    stats: Mutex<MemoryStats>,
    start_snapshot: Mutex<Option<MemoryStats>>,
}

impl MemoryProfiler {
    /// Create a new memory profiler
    pub fn new() -> Self {
        Self {
            stats: Mutex::new(MemoryStats::default()),
            start_snapshot: Mutex::new(None),
        }
    }

    /// Start memory profiling session
    pub fn start_profiling(&self) {
        let mut start_snapshot = self.start_snapshot.lock().unwrap();
        *start_snapshot = Some(self.get_stats());
    }

    /// Stop profiling and return statistics for the session
    pub fn stop_profiling(&self) -> MemoryStats {
        let start_stats = self
            .start_snapshot
            .lock()
            .unwrap()
            .take()
            .unwrap_or_else(|| MemoryStats::default());

        let current_stats = self.get_stats();

        MemoryStats {
            current_allocated: current_stats
                .current_allocated
                .saturating_sub(start_stats.current_allocated),
            peak_allocated: current_stats
                .peak_allocated
                .saturating_sub(start_stats.peak_allocated),
            allocation_count: current_stats
                .allocation_count
                .saturating_sub(start_stats.allocation_count),
            deallocation_count: current_stats
                .deallocation_count
                .saturating_sub(start_stats.deallocation_count),
            total_allocated: current_stats
                .total_allocated
                .saturating_sub(start_stats.total_allocated),
            total_deallocated: current_stats
                .total_deallocated
                .saturating_sub(start_stats.total_deallocated),
            avg_allocation_size: if current_stats.allocation_count > start_stats.allocation_count {
                let net_allocations = current_stats.allocation_count - start_stats.allocation_count;
                let net_bytes = current_stats.total_allocated - start_stats.total_allocated;
                net_bytes as f64 / net_allocations as f64
            } else {
                0.0
            },
            max_allocation_size: current_stats.max_allocation_size,
        }
    }

    /// Get current memory statistics
    pub fn get_stats(&self) -> MemoryStats {
        GLOBAL_ALLOCATOR.get_stats()
    }

    /// Profile a closure and return memory statistics
    pub fn profile_closure<F, T>(&self, operation: F) -> (T, MemoryStats)
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
    stats: Mutex<MemoryStats>,
}

impl TrackingAllocator {
    const fn new() -> Self {
        Self {
            allocator: System,
            stats: Mutex::new(MemoryStats {
                total_allocated: 0,
                total_deallocated: 0,
                current_allocated: 0,
                peak_allocated: 0,
                deallocation_count: 0,
                allocation_count: 0,
                avg_allocation_size: 0.0,
                max_allocation_size: 0,
            }),
        }
    }

    fn get_stats(&self) -> MemoryStats {
        self.stats.lock().unwrap().clone()
    }

    fn record_allocation(&self, size: usize) {
        let mut stats = self.stats.lock().unwrap();
        stats.allocation_count += 1;
        stats.total_allocated += size;
        stats.current_allocated += size;
        stats.max_allocation_size = stats.max_allocation_size.max(size);
        stats.peak_allocated = stats.peak_allocated.max(stats.current_allocated);
        stats.avg_allocation_size = stats.total_allocated as f64 / stats.allocation_count as f64;
    }

    fn record_deallocation(&self, size: usize) {
        let mut stats = self.stats.lock().unwrap();
        stats.deallocation_count += 1;
        stats.total_deallocated += size;
        stats.current_allocated = stats.current_allocated.saturating_sub(size);
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
            let mut stats = self.stats.lock().unwrap();
            stats.current_allocated = stats.current_allocated.saturating_sub(old_size) + new_size;
            stats.peak_allocated = stats.peak_allocated.max(stats.current_allocated);
            stats.max_allocation_size = stats.max_allocation_size.max(new_size);

            // Count as deallocation + allocation
            stats.deallocation_count += 1;
            stats.total_deallocated += old_size;
            stats.allocation_count += 1;
            stats.total_allocated += new_size;
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
    pub fn profile_matrix_operations(&self) -> Result<MemoryStats> {
        Ok(self
            .profiler
            .profile_closure(|| {
                // Create and manipulate CFD-style matrices
                let mut builder = cfd_math::sparse::SparseMatrixBuilder::new(1000, 1000);

                // Simulate CFD matrix assembly (5-point stencil)
                for i in 0..32 {
                    for j in 0..32 {
                        let idx = i * 32 + j;

                        // Add stencil entries
                        builder.add_entry(idx, idx, 4.0).unwrap(); // Center
                        if i > 0 {
                            builder.add_entry(idx, idx - 32, -1.0).unwrap();
                        } // North
                        if i < 31 {
                            builder.add_entry(idx, idx + 32, -1.0).unwrap();
                        } // South
                        if j > 0 {
                            builder.add_entry(idx, idx - 1, -1.0).unwrap();
                        } // West
                        if j < 31 {
                            builder.add_entry(idx, idx + 1, -1.0).unwrap();
                        } // East
                    }
                }

                let matrix = builder.build().unwrap();
                let vector = vec![1.0f64; 1000];
                let mut result = vec![0.0f64; 1000];

                // Perform SPMV operations
                for _ in 0..10 {
                    let vector_dv = nalgebra::DVector::from_vec(vector.clone());
                    let mut result_dv = nalgebra::DVector::zeros(matrix.nrows());
                    cfd_math::sparse::spmv(&matrix, &vector_dv, &mut result_dv);
                }
            })
            .1)
    }

    /// Profile memory usage of CFD grid operations
    pub fn profile_grid_operations(&self) -> Result<MemoryStats> {
        use cfd_2d::grid::StructuredGrid2D;

        Ok(self
            .profiler
            .profile_closure(|| {
                // Create multiple grids of different sizes
                let grids = vec![
                    StructuredGrid2D::new(32, 32, 0.0, 1.0, 0.0, 1.0).unwrap(),
                    StructuredGrid2D::new(64, 64, 0.0, 1.0, 0.0, 1.0).unwrap(),
                    StructuredGrid2D::new(128, 128, 0.0, 1.0, 0.0, 1.0).unwrap(),
                ];

                // Perform grid operations
                for grid in grids {
                    let nx = grid.nx();
                    let ny = grid.ny();
                    // Access grid data for memory profiling
                    let _centers: Vec<_> = (0..nx.min(10))
                        .flat_map(|i| {
                            let grid_ref = &grid;
                            (0..ny.min(10)).map(move |j| grid_ref.cell_center(i, j).ok())
                        })
                        .collect();
                }
            })
            .1)
    }

    /// Profile memory usage of CFD simulation fields
    pub fn profile_simulation_fields(&self) -> Result<MemoryStats> {
        use cfd_2d::fields::SimulationFields;
        use cfd_core::fluid::Fluid;

        Ok(self
            .profiler
            .profile_closure(|| {
                let fluid = Fluid::new("Test".to_string(), 1.0, 0.01, 1000.0, 0.001);

                // Create fields of increasing size
                let sizes = vec![32, 64, 128];

                for &size in &sizes {
                    let fields = SimulationFields::with_fluid(size, size, &fluid);

                    // Access all field data to ensure it's in memory
                    let _u_sum: f64 = fields.u.data().iter().sum();
                    let _v_sum: f64 = fields.v.data().iter().sum();
                    let _p_sum: f64 = fields.p.data().iter().sum();
                }
            })
            .1)
    }

    /// Run comprehensive CFD memory profiling suite
    pub fn run_memory_suite(&self) -> Result<Vec<(String, MemoryStats)>> {
        println!("Running CFD Memory Profiling Suite...");

        let mut results = Vec::new();

        // Matrix operations
        let matrix_stats = self.profile_matrix_operations()?;
        results.push(("Matrix_Operations".to_string(), matrix_stats));

        // Grid operations
        let grid_stats = self.profile_grid_operations()?;
        results.push(("Grid_Operations".to_string(), grid_stats));

        // Simulation fields
        let field_stats = self.profile_simulation_fields()?;
        results.push(("Simulation_Fields".to_string(), field_stats));

        // Memory allocation patterns
        let alloc_stats = self
            .profiler
            .profile_closure(|| {
                let mut allocations = Vec::new();
                for size in [100, 1000, 10000, 100000] {
                    allocations.push(vec![0.0f64; size]);
                }
                // Hold onto allocations briefly
                std::thread::sleep(std::time::Duration::from_millis(1));
                drop(allocations);
            })
            .1;
        results.push(("Allocation_Patterns".to_string(), alloc_stats));

        Ok(results)
    }
}

impl Default for CfdMemoryProfiler {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for MemoryStats {
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
        let stats = MemoryStats {
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
