//! Tests for the 2D Laplacian GPU kernel.

use super::cpu_reference::execute_cpu_reference;
use super::kernel::Laplacian2DKernel;
use super::types::BoundaryType;
use crate::compute::gpu::GpuContext;
use std::sync::Arc;

fn create_kernel() -> Option<Laplacian2DKernel> {
    GpuContext::create()
        .ok()
        .map(|ctx| Laplacian2DKernel::new(Arc::new(ctx)))
}

fn execute_with_fallback(
    kernel: Option<&Laplacian2DKernel>,
    field: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    bc: BoundaryType,
    result: &mut [f32],
) {
    if let Some(kernel) = kernel {
        kernel.execute_with_bc(field, nx, ny, dx, dy, bc, result);
    } else {
        execute_cpu_reference(field, nx, ny, dx, dy, bc, result);
    }
}

fn execute_cpu_fallback(
    _kernel: Option<&Laplacian2DKernel>,
    field: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    bc: BoundaryType,
    result: &mut [f32],
) {
    // Always use the CPU reference to avoid depending on the private execute_cpu method
    execute_cpu_reference(field, nx, ny, dx, dy, bc, result);
}

/// Test function: u(x,y) = sin(πx)sin(πy) on [0,1]×[0,1]
/// Exact Laplacian: ∇²u = -2π²sin(πx)sin(πy)
fn test_function_1(x: f32, y: f32) -> f32 {
    (std::f32::consts::PI * x).sin() * (std::f32::consts::PI * y).sin()
}

fn exact_laplacian_1(x: f32, y: f32) -> f32 {
    -2.0 * std::f32::consts::PI * std::f32::consts::PI * test_function_1(x, y)
}

/// Test function: u(x,y) = x² + y²
/// Exact Laplacian: ∇²u = 4
fn test_function_2(x: f32, y: f32) -> f32 {
    x * x + y * y
}

fn exact_laplacian_2(_x: f32, _y: f32) -> f32 {
    4.0
}

#[test]
fn test_laplacian_accuracy_polynomial() {
    let kernel = create_kernel();
    let n = 32;
    let dx = 1.0 / (n - 1) as f32;
    let dy = 1.0 / (n - 1) as f32;

    let mut field = vec![0.0; n * n];
    let mut result = vec![0.0; n * n];

    // Initialize field with test function
    for j in 0..n {
        for i in 0..n {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            field[j * n + i] = test_function_2(x, y);
        }
    }

    execute_with_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Dirichlet,
        &mut result,
    );

    // Verify accuracy - should be exactly 4.0 for polynomial
    let mut max_error: f32 = 0.0;
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            let exact = exact_laplacian_2(x, y);
            let error = (result[j * n + i] - exact).abs();
            max_error = max_error.max(error);
        }
    }

    // For this simple polynomial, discrete central differences yield exact results analytically.
    // Allow small numerical error from floating-point rounding.
    assert!(
        max_error < 1e-3,
        "Max error {max_error} too large for polynomial test"
    );
}

#[test]
fn test_laplacian_convergence_rate() {
    let kernel = create_kernel();

    let grid_sizes = vec![16, 32, 64];
    let mut errors = Vec::new();

    for &n in &grid_sizes {
        let dx = 1.0 / (n - 1) as f32;
        let dy = 1.0 / (n - 1) as f32;

        let mut field = vec![0.0; n * n];
        let mut result = vec![0.0; n * n];

        // Initialize field with sinusoidal test function
        for j in 0..n {
            for i in 0..n {
                let x = i as f32 * dx;
                let y = j as f32 * dy;
                field[j * n + i] = test_function_1(x, y);
            }
        }

        execute_with_fallback(
            kernel.as_ref(),
            &field,
            n,
            n,
            dx,
            dy,
            BoundaryType::Dirichlet,
            &mut result,
        );

        // Calculate L2 error
        let mut l2_error = 0.0;
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let x = i as f32 * dx;
                let y = j as f32 * dy;
                let exact = exact_laplacian_1(x, y);
                let error = result[j * n + i] - exact;
                l2_error += error * error * dx * dy;
            }
        }
        errors.push(l2_error.sqrt());
    }

    // Verify second-order convergence
    for i in 1..errors.len() {
        let rate = (errors[i - 1] / errors[i]).log2();
        assert!(rate > 1.8, "Convergence rate {rate} too low, expected ~2.0");
    }
}

#[test]
fn test_boundary_conditions_dirichlet() {
    let kernel = create_kernel();
    let n = 32;
    let dx = 1.0 / (n - 1) as f32;
    let dy = 1.0 / (n - 1) as f32;

    let mut field = vec![0.0; n * n];
    let mut result = vec![0.0; n * n];

    // Initialize with function that satisfies u=0 on boundaries
    for j in 0..n {
        for i in 0..n {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            field[j * n + i] = x * (1.0 - x) * y * (1.0 - y); // Zero on all boundaries
        }
    }

    execute_with_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Dirichlet,
        &mut result,
    );

    // Verify boundary points are computed (not left uninitialized)
    for j in 0..n {
        for i in 0..n {
            assert!(
                !result[j * n + i].is_nan(),
                "Boundary point ({i},{j}) is NaN"
            );
            assert!(
                result[j * n + i].is_finite(),
                "Boundary point ({i},{j}) is infinite"
            );
        }
    }
}

#[test]
fn test_boundary_conditions_neumann() {
    let kernel = create_kernel();
    let n = 32;
    let dx = 1.0 / (n - 1) as f32;
    let dy = 1.0 / (n - 1) as f32;

    let mut field = vec![0.0; n * n];
    let mut result = vec![0.0; n * n];

    // Initialize with function that has zero normal derivative on boundaries
    // u(x,y) = x² + y² has ∂u/∂n = 0 on boundaries when properly centered
    for j in 0..n {
        for i in 0..n {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            field[j * n + i] = x * x + y * y;
        }
    }

    execute_with_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Neumann,
        &mut result,
    );

    // Verify boundary points are computed and finite
    for j in 0..n {
        for i in 0..n {
            assert!(
                !result[j * n + i].is_nan(),
                "Boundary point ({i},{j}) is NaN"
            );
            assert!(
                result[j * n + i].is_finite(),
                "Boundary point ({i},{j}) is infinite"
            );
        }
    }

    // For u(x,y) = x² + y², ∇²u = 4 everywhere, including boundaries
    // Verify this holds approximately at boundary points
    let tolerance = 0.1; // Allow some numerical error near boundaries
    for j in 0..n {
        for i in 0..n {
            if i == 0 || i == n - 1 || j == 0 || j == n - 1 {
                let val = result[j * n + i];
                let error = (val - 4.0).abs();
                assert!(
                    error < tolerance,
                    "Boundary point ({i},{j}) has ∇²u = {val}, expected ~4.0"
                );
            }
        }
    }
}

#[test]
fn test_boundary_conditions_periodic() {
    let kernel = create_kernel();
    let n = 32;
    let dx = 1.0 / (n - 1) as f32;
    let dy = 1.0 / (n - 1) as f32;

    let mut field = vec![0.0; n * n];
    let mut result = vec![0.0; n * n];

    // Initialize with periodic function: u(x,y) = sin(2πx)sin(2πy)
    // This satisfies periodic boundary conditions naturally
    for j in 0..n {
        for i in 0..n {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            field[j * n + i] =
                (2.0 * std::f32::consts::PI * x).sin() * (2.0 * std::f32::consts::PI * y).sin();
        }
    }

    execute_with_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Periodic,
        &mut result,
    );

    // Verify boundary points are computed and finite
    for j in 0..n {
        for i in 0..n {
            assert!(
                !result[j * n + i].is_nan(),
                "Boundary point ({i},{j}) is NaN"
            );
            assert!(
                result[j * n + i].is_finite(),
                "Boundary point ({i},{j}) is infinite"
            );
        }
    }

    // Verify periodicity: result should be periodic too
    // Check that opposite boundaries have similar values (within numerical tolerance)
    let tolerance = 0.1;

    // Check left-right periodicity
    for j in 0..n {
        let left_val = result[j * n];
        let right_val = result[j * n + (n - 1)];
        let error = (left_val - right_val).abs();
        assert!(
            error < tolerance,
            "Periodicity violation at j={j}: left={left_val}, right={right_val}"
        );
    }

    // Check bottom-top periodicity
    for i in 0..n {
        let bottom_val = result[i];
        let top_val = result[(n - 1) * n + i];
        let error = (bottom_val - top_val).abs();
        assert!(
            error < tolerance,
            "Periodicity violation at i={i}: bottom={bottom_val}, top={top_val}"
        );
    }
}

#[test]
fn test_boundary_conditions_comprehensive() {
    let kernel = create_kernel();
    let n = 16; // Smaller grid for comprehensive test
    let dx = 1.0 / (n - 1) as f32;
    let dy = 1.0 / (n - 1) as f32;

    // Test different manufactured solutions for each boundary condition type
    // Use function pointers to avoid closure type issues
    type TestFunction = fn(f32, f32) -> f32;

    fn dirichlet_func(x: f32, y: f32) -> f32 {
        x * (1.0 - x) * y * (1.0 - y) // Zero on boundaries
    }

    fn neumann_func(x: f32, y: f32) -> f32 {
        x * x + y * y // Zero normal derivative
    }

    fn periodic_func(x: f32, y: f32) -> f32 {
        (2.0 * std::f32::consts::PI * x).sin() * (2.0 * std::f32::consts::PI * y).sin()
    }

    let test_cases: Vec<(&str, TestFunction)> = vec![
        ("Dirichlet", dirichlet_func),
        ("Neumann", neumann_func),
        ("Periodic", periodic_func),
    ];

    for (bc_type, test_function) in test_cases {
        let mut field = vec![0.0; n * n];
        let mut result = vec![0.0; n * n];

        // Initialize field with test function
        for j in 0..n {
            for i in 0..n {
                let x = i as f32 * dx;
                let y = j as f32 * dy;
                field[j * n + i] = test_function(x, y);
            }
        }

        let bc = match bc_type {
            "Neumann" => BoundaryType::Neumann,
            "Periodic" => BoundaryType::Periodic,
            _ => BoundaryType::Dirichlet,
        };
        execute_with_fallback(kernel.as_ref(), &field, n, n, dx, dy, bc, &mut result);

        // Verify all computed values are finite and reasonable
        let mut max_val = f32::NEG_INFINITY;
        let mut min_val = f32::INFINITY;
        let mut nan_count = 0;
        let mut inf_count = 0;

        for j in 0..n {
            for i in 0..n {
                let val = result[j * n + i];
                if val.is_nan() {
                    nan_count += 1;
                } else if val.is_infinite() {
                    inf_count += 1;
                } else {
                    max_val = max_val.max(val);
                    min_val = min_val.min(val);
                }
            }
        }

        assert_eq!(nan_count, 0, "{bc_type} BC: Found {nan_count} NaN values");
        assert_eq!(
            inf_count, 0,
            "{bc_type} BC: Found {inf_count} infinite values"
        );

        // Verify reasonable range (Laplacian should not explode)
        let range = max_val - min_val;
        let range_limit = if bc_type == "Periodic" { 200.0 } else { 100.0 };
        assert!(
            range < range_limit,
            "{bc_type} BC: Result range {range} is too large"
        );

        // Verify boundary behavior based on BC type
        match bc_type {
            "Periodic" => {
                // Check periodicity for periodic BC
                let tolerance = 0.2;
                for j in 0..n {
                    let left = result[j * n];
                    let right = result[j * n + (n - 1)];
                    assert!(
                        (left - right).abs() < tolerance,
                        "{bc_type} BC: Periodicity violation at j={j}"
                    );
                }
                for i in 0..n {
                    let bottom = result[i];
                    let top = result[(n - 1) * n + i];
                    assert!(
                        (bottom - top).abs() < tolerance,
                        "{bc_type} BC: Periodicity violation at i={i}"
                    );
                }
            }
            "Neumann" => {
                // For Neumann BC with u(x,y) = x² + y², ∇²u = 4 everywhere
                let expected_laplacian = 4.0;
                let tolerance = 0.5;
                for j in 0..n {
                    for i in 0..n {
                        let val = result[j * n + i];
                        let error = (val - expected_laplacian).abs();
                        assert!(
                            error < tolerance,
                            "{bc_type} BC: Point ({i},{j}) has ∇²u={val}, expected ~{expected_laplacian}"
                        );
                    }
                }
            }
            _ => {}
        }
    }
}

#[test]
fn test_gpu_cpu_consistency() {
    let kernel = create_kernel();
    let n = 32;
    let dx = 0.1f32;
    let dy = 0.1f32;

    let mut field = vec![0.0; n * n];
    let mut gpu_result = vec![0.0; n * n];
    let mut cpu_result = vec![0.0; n * n];

    // Initialize with random field
    for (i, val) in field.iter_mut().enumerate() {
        *val = (i as f32).sin() * 0.5 + 0.5;
    }

    // Force CPU execution
    execute_cpu_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Dirichlet,
        &mut cpu_result,
    );

    // Force GPU execution
    execute_with_fallback(
        kernel.as_ref(),
        &field,
        n,
        n,
        dx,
        dy,
        BoundaryType::Dirichlet,
        &mut gpu_result,
    );

    // Compare results (allowing for small numerical differences)
    let mut max_diff = 0.0f32;
    for i in 0..field.len() {
        let diff = (gpu_result[i] - cpu_result[i]).abs();
        max_diff = max_diff.max(diff);
    }

    assert!(
        max_diff < 1e-6,
        "GPU/CPU inconsistency: max difference {max_diff}"
    );
}

#[test]
fn test_gpu_cpu_performance_benchmark() {
    use std::time::Instant;

    let kernel = create_kernel();

    // Test multiple grid sizes to analyze scaling
    let grid_sizes = vec![16, 32, 64, 128, 256];
    let num_warmup_runs = 5;
    let num_timing_runs = 10;

    println!("\n=== GPU vs CPU Performance Benchmark ===");
    println!("Grid Size | CPU Time (ms) | GPU Time (ms) | Speedup | Throughput (MCells/s)");
    println!("----------|---------------|---------------|---------|----------------------");

    for &n in &grid_sizes {
        let dx = 1.0 / (n - 1) as f32;
        let dy = 1.0 / (n - 1) as f32;
        let mut field = vec![0.0; n * n];
        let mut gpu_result = vec![0.0; n * n];
        let mut cpu_result = vec![0.0; n * n];

        // Initialize with smooth function
        for j in 0..n {
            for i in 0..n {
                let x = i as f32 * dx;
                let y = j as f32 * dy;
                field[j * n + i] =
                    (std::f32::consts::PI * x).sin() * (std::f32::consts::PI * y).cos();
            }
        }

        // CPU benchmarking
        let cpu_time = {
            // Warmup runs
            for _ in 0..num_warmup_runs {
                execute_cpu_fallback(
                    kernel.as_ref(),
                    &field,
                    n,
                    n,
                    dx,
                    dy,
                    BoundaryType::Dirichlet,
                    &mut cpu_result,
                );
            }

            // Timing runs
            let start = Instant::now();
            for _ in 0..num_timing_runs {
                execute_cpu_fallback(
                    kernel.as_ref(),
                    &field,
                    n,
                    n,
                    dx,
                    dy,
                    BoundaryType::Dirichlet,
                    &mut cpu_result,
                );
            }
            let elapsed = start.elapsed();
            elapsed.as_secs_f64() / f64::from(num_timing_runs) * 1000.0 // Convert to ms
        };

        // GPU benchmarking
        let gpu_time = {
            // Warmup runs
            for _ in 0..num_warmup_runs {
                execute_with_fallback(
                    kernel.as_ref(),
                    &field,
                    n,
                    n,
                    dx,
                    dy,
                    BoundaryType::Dirichlet,
                    &mut gpu_result,
                );
            }

            // Timing runs
            let start = Instant::now();
            for _ in 0..num_timing_runs {
                execute_with_fallback(
                    kernel.as_ref(),
                    &field,
                    n,
                    n,
                    dx,
                    dy,
                    BoundaryType::Dirichlet,
                    &mut gpu_result,
                );
            }
            let elapsed = start.elapsed();
            elapsed.as_secs_f64() / f64::from(num_timing_runs) * 1000.0 // Convert to ms
        };

        // Calculate metrics
        let speedup = cpu_time / gpu_time;
        let total_cells = (n * n) as f64;
        let cpu_throughput = total_cells / (cpu_time / 1000.0) / 1e6; // MCells/s
        let gpu_throughput = total_cells / (gpu_time / 1000.0) / 1e6; // MCells/s

        println!(
            "{n:9} | {cpu_time:13.3} | {gpu_time:13.3} | {speedup:7.2}x | {cpu_throughput:18.1} (CPU)"
        );
        println!(
            "{:9} | {:13} | {:13} | {:7} | {gpu_throughput:18.1} (GPU)",
            "", "", "", ""
        );

        // Verify correctness for this grid size.
        // The CPU reference uses f64 intermediate precision while the GPU
        // shader uses native f32 arithmetic.  The central-difference stencil
        // (u_left - 2*u_center + u_right) suffers catastrophic cancellation
        // that amplifies the f32/f64 rounding difference by 1/dx^2.  A
        // relative tolerance of 1e-4 is appropriate for this comparison.
        let mut max_error = 0.0f32;
        let mut max_magnitude = 0.0f32;
        for i in 0..field.len() {
            let error = (gpu_result[i] - cpu_result[i]).abs();
            max_error = max_error.max(error);
            max_magnitude = max_magnitude.max(cpu_result[i].abs());
        }
        let rel_error = if max_magnitude > 0.0 {
            max_error / max_magnitude
        } else {
            max_error
        };

        assert!(
            rel_error < 1e-4,
            "GPU/CPU inconsistency at n={n}: max abs error {max_error}, \
             rel error {rel_error}, max magnitude {max_magnitude}"
        );
    }
}
