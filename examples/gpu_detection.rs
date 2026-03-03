//! GPU Detection Example
//!
//! This example demonstrates automatic GPU detection with support for:
//! - Discrete GPUs (NVIDIA, AMD)
//! - Integrated GPUs (Intel, AMD APU)
//! - Software fallback (SwiftShader, LLVMpipe)

use cfd_core::compute::backend::ComputeCapability;
use cfd_core::compute::dispatch::ComputeDispatcher;
use cfd_core::compute::traits::ComputeBackend;

fn main() {
    println!("=== GPU/SIMD Detection Test ===\n");

    // Detect system capabilities
    let capabilities = ComputeCapability::detect();

    println!("Available backends: {:?}", capabilities.backends);
    println!("Preferred backend: {:?}", capabilities.preferred_backend);
    println!("Compute units: {}", capabilities.compute_units);
    println!("Available memory: {} MB", capabilities.available_memory / (1024 * 1024));

    // Initialize unified compute dispatcher
    match ComputeDispatcher::new() {
        Ok(dispatcher) => {
            println!("\nCompute initialization successful!");

            match dispatcher.current_backend() {
                ComputeBackend::Gpu => {
                    println!("GPU acceleration is active");
                    println!("   Your system's GPU (discrete or integrated) will be used");
                }
                ComputeBackend::Simd => {
                    println!("SIMD acceleration is active");
                    println!("   Hardware vector instructions will be used");
                }
                ComputeBackend::Cpu => {
                    println!("CPU fallback is active");
                    println!("   Scalar CPU operations will be used");
                }
                ComputeBackend::Hybrid => {
                    println!("Hybrid (GPU+CPU) acceleration is active");
                }
            }

            // Run a simple test using raw vectors (no kernel dispatch)
            println!("\n=== Performance Test ===");
            let size = 1_000_000;
            let a = vec![1.0f32; size];
            let b = vec![2.0f32; size];
            let mut result = vec![0.0f32; size];

            let start = std::time::Instant::now();
            for i in 0..size {
                result[i] = a[i] + b[i];
            }
            let duration = start.elapsed();

            println!("Vector addition (1M elements): {:?}", duration);

            // Verify correctness
            let sample_check = result.iter().take(10).all(|&x| (x - 3.0).abs() < 1e-6);
            if sample_check {
                println!("Computation correct");
            } else {
                println!("Computation error");
            }

            // Estimate throughput
            let gb_processed = (size * 3 * 4) as f64 / 1e9; // 3 arrays, 4 bytes each
            let throughput = gb_processed / duration.as_secs_f64();
            println!("Throughput: {:.2} GB/s", throughput);

            println!("\n=== System Capabilities ===");
            #[cfg(feature = "gpu")]
            {
                println!("GPU support: Enabled (supports integrated graphics)");
            }
            #[cfg(not(feature = "gpu"))]
            {
                println!("GPU support: Disabled");
            }

            // Check SIMD capabilities
            use cfd_math::simd::SimdCapability;
            let simd_cap = SimdCapability::detect();
            println!("SIMD capability: {:?}", simd_cap);

            // Report backend selection for various problem sizes
            println!("\n=== Backend Selection by Problem Size ===");
            for &problem_size in &[100, 10_000, 1_000_000] {
                let selected = dispatcher.capabilities().select_backend(problem_size);
                println!("  Problem size {:>10}: {:?}", problem_size, selected);
            }
        }
        Err(e) => {
            println!("Failed to initialize compute: {}", e);
        }
    }
}
