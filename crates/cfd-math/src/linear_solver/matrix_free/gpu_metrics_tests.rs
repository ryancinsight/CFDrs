#![cfg(feature = "gpu")]

use super::gpu_compute::{GpuComputeContext};
use super::gpu_operators::{GpuLaplacianOperator2D, BoundaryType};

#[tokio::test]
async fn test_gpu_dispatch_metrics_present_or_skip() {
    // Try to create a GPU context; if unavailable, skip
    let ctx = match GpuComputeContext::new().await {
        Ok(c) => c,
        Err(_) => return, // Skip on systems without GPU adapter
    };

    let nx = 32usize;
    let ny = 32usize;
    let dx = 1.0f32 / (nx as f32 - 1.0);
    let dy = 1.0f32 / (ny as f32 - 1.0);

    let op = GpuLaplacianOperator2D::new(
        std::sync::Arc::new(ctx),
        nx,
        ny,
        dx,
        dy,
        BoundaryType::Dirichlet,
    );

    let mut field = vec![0.0f32; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let x = i as f32 * dx;
            let y = j as f32 * dy;
            field[j * nx + i] = x * x + y * y;
        }
    }
    let mut lap = vec![0.0f32; nx * ny];
    let metrics = op.apply_gpu_with_metrics(&field, &mut lap).await.unwrap();
    assert!(metrics.duration_ms >= 0.0);
}
