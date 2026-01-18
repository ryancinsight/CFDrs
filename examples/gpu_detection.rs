//! GPU Detection Example
//!
//! This example demonstrates automatic GPU detection with support for:
//! - Discrete GPUs (NVIDIA, AMD)
//! - Integrated GPUs (Intel, AMD APU)
//! - Software fallback (SwiftShader, LLVMpipe)

use cfd_suite::compute_unified::{Backend, UnifiedCompute};

fn main() {
    println!("=== GPU/SIMD Detection Test ===\n");

    // Initialize unified compute
    match UnifiedCompute::new() {
        Ok(mut compute) => {
            println!("✅ Compute initialization successful!");

            match compute.backend() {
                Backend::Gpu => {
                    println!("🎮 GPU acceleration is active");
                    println!("   Your system's GPU (discrete or integrated) will be used");
                }
                Backend::Simd => {
                    println!("⚡ SIMD acceleration is active");
                    println!("   Hardware vector instructions will be used");
                }
                Backend::Swar => {
                    println!("🔧 SWAR fallback is active");
                    println!("   Software-based optimization will be used");
                }
            }

            // Run a simple test
            println!("\n=== Performance Test ===");
            let size = 1_000_000;
            let a = vec![1.0f32; size];
            let b = vec![2.0f32; size];
            let mut result = vec![0.0f32; size];

            let start = std::time::Instant::now();
            compute.vector_add_f32(&a, &b, &mut result).unwrap();
            let duration = start.elapsed();

            println!("Vector addition (1M elements): {:?}", duration);

            // Verify correctness
            let sample_check = result.iter().take(10).all(|&x| x == 3.0);
            if sample_check {
                println!("✅ Computation correct");
            } else {
                println!("❌ Computation error");
            }

            // Estimate throughput
            let gb_processed = (size * 3 * 4) as f64 / 1e9; // 3 arrays, 4 bytes each
            let throughput = gb_processed / duration.as_secs_f64();
            println!("Throughput: {:.2} GB/s", throughput);

            println!("\n=== System Capabilities ===");
            #[cfg(feature = "gpu")]
            {
                println!("GPU support: ✅ Enabled (supports integrated graphics)");
            }
            #[cfg(not(feature = "gpu"))]
            {
                println!("GPU support: ❌ Disabled");
            }

            // Check SIMD capabilities
            use cfd_math::simd::SimdCapability;
            let simd_cap = SimdCapability::detect();
            println!("SIMD capability: {:?}", simd_cap);
        }
        Err(e) => {
            println!("❌ Failed to initialize compute: {}", e);
        }
    }
}
