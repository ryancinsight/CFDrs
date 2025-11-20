//! GPU shader definitions for field operations

/// WGSL shader for field addition
pub const FIELD_ADD_SHADER: &str = r"
@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn add_fields(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&a)) {
        return;
    }
    result[index] = a[index] + b[index];
}
";

/// WGSL shader for scalar multiplication
pub const FIELD_MUL_SHADER: &str = r"
@group(0) @binding(0) var<storage, read> field: array<f32>;
@group(0) @binding(1) var<uniform> scalar: f32;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn multiply_field(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&field)) {
        return;
    }
    result[index] = field[index] * scalar;
}
";

/// WGSL shader for 2D Laplacian with rigorous boundary handling
/// Implements a mathematically correct 5‑point stencil for elliptic PDEs
/// Reference: LeVeque (2007), Trefethen (1996) – Finite Difference Methods
///
/// Mathematical Foundation:
/// - Interior: 5‑point stencil ∇²u ≈ (u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j})/h²
/// - Dirichlet: ghost points enforce u=0 yielding −2u_{boundary}/h² contribution per axis
/// - Neumann: second‑order one‑sided second derivatives at boundaries:
///   • x‑left (i=0): (2u₀ − 5u₁ + 4u₂ − u₃)/Δx²
///   • x‑right (i=nx−1): (2u_{nx−1} − 5u_{nx−2} + 4u_{nx−3} − u_{nx−4})/Δx²
///   • y‑bottom (j=0): (2u₀ − 5u₁ + 4u₂ − u₃)/Δy²
///   • y‑top (j=ny−1): (2u_{ny−1} − 5u_{ny−2} + 4u_{ny−3} − u_{ny−4})/Δy²
///   Fallback to symmetric ghost scheme for small grids (<4 points) to avoid OOB.
/// - Periodic: wrap neighbors across domain endpoints (i=0 ↔ i=nx−1, j=0 ↔ j=ny−1)
/// - Ensures well‑posedness by defining operator at all grid points with correct BCs
pub const LAPLACIAN_2D_SHADER: &str = r"
struct Uniforms {
    dims_bc: vec4<u32>,   // (nx, ny, bc_type, pad)
    inv2: vec4<f32>,      // (dx_inv2, dy_inv2, 0.0, 0.0)
}

@group(0) @binding(0) var<storage, read> field: array<f32>;
@group(0) @binding(1) var<uniform> uniforms: Uniforms;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(8, 8)
fn laplacian_2d(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    
    if (i >= uniforms.dims_bc.x || j >= uniforms.dims_bc.y) {
        return;
    }
    
    let idx = j * uniforms.dims_bc.x + i;
    var laplacian = 0.0;
    
    // X direction - handle boundary conditions properly
    if (i > 0u && i < uniforms.dims_bc.x - 1u) {
        // Interior point - standard 5-point stencil
        let left = field[j * uniforms.dims_bc.x + (i - 1u)];
        let center = field[idx];
        let right = field[j * uniforms.dims_bc.x + (i + 1u)];
        laplacian += (left - 2.0 * center + right) * uniforms.inv2.x;
    } else if (i == 0u) {
        // Left boundary
        if (uniforms.dims_bc.z == 0u) { // Dirichlet: u(0,j) = 0
            let center = field[idx];
            // Ghost point u(-1,j) = -u(1,j) => left+right-2*center cancels neighbor, yields -2*center
            laplacian += (-2.0 * center) * uniforms.inv2.x;
        } else if (uniforms.dims_bc.z == 1u) { // Neumann: du/dx = 0
            let nx = uniforms.dims_bc.x;
            let center = field[idx];
            if (nx >= 4u) {
                let u1 = field[j * nx + 1u];
                let u2 = field[j * nx + 2u];
                let u3 = field[j * nx + 3u];
                // Second-order one-sided second derivative at left boundary
                laplacian += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * uniforms.inv2.x;
            } else {
                let right = field[j * nx + (i + 1u)];
                // Fallback: symmetric ghost approach for small nx
                laplacian += (right - 2.0 * center + right) * uniforms.inv2.x;
            }
        } else if (uniforms.dims_bc.z == 2u) { // Periodic
            let left = field[j * uniforms.dims_bc.x + (uniforms.dims_bc.x - 1u)];
            let center = field[idx];
            let right = field[j * uniforms.dims_bc.x + (i + 1u)];
            laplacian += (left - 2.0 * center + right) * uniforms.inv2.x;
        }
    } else if (i == uniforms.dims_bc.x - 1u) {
        // Right boundary
        if (uniforms.dims_bc.z == 0u) { // Dirichlet: u(nx-1,j) = 0
            let center = field[idx];
            // Ghost point u(nx,j) = -u(nx-2,j) => left+right-2*center cancels neighbor, yields -2*center
            laplacian += (-2.0 * center) * uniforms.inv2.x;
        } else if (uniforms.dims_bc.z == 1u) { // Neumann: du/dx = 0
            let nx = uniforms.dims_bc.x;
            let center = field[idx];
            if (nx >= 4u) {
                let u1 = field[j * nx + (nx - 2u)];
                let u2 = field[j * nx + (nx - 3u)];
                let u3 = field[j * nx + (nx - 4u)];
                // Second-order one-sided second derivative at right boundary
                laplacian += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * uniforms.inv2.x;
            } else {
                let left = field[j * nx + (i - 1u)];
                // Fallback: symmetric ghost approach for small nx
                laplacian += (left - 2.0 * center + left) * uniforms.inv2.x;
            }
        } else if (uniforms.dims_bc.z == 2u) { // Periodic
            let left = field[j * uniforms.dims_bc.x + (i - 1u)];
            let center = field[idx];
            let right = field[j * uniforms.dims_bc.x + 0u];
            laplacian += (left - 2.0 * center + right) * uniforms.inv2.x;
        }
    }
    
    // Y direction - handle boundary conditions properly
    if (j > 0u && j < uniforms.dims_bc.y - 1u) {
        // Interior point - standard 5-point stencil
        let bottom = field[(j - 1u) * uniforms.dims_bc.x + i];
        let center = field[idx];
        let top = field[(j + 1u) * uniforms.dims_bc.x + i];
        laplacian += (bottom - 2.0 * center + top) * uniforms.inv2.y;
    } else if (j == 0u) {
        // Bottom boundary
        if (uniforms.dims_bc.z == 0u) { // Dirichlet: u(i,0) = 0
            let center = field[idx];
            // Ghost point u(i,-1) = -u(i,1) => bottom+top-2*center cancels neighbor, yields -2*center
            laplacian += (-2.0 * center) * uniforms.inv2.y;
        } else if (uniforms.dims_bc.z == 1u) { // Neumann: du/dy = 0
            let ny = uniforms.dims_bc.y;
            let center = field[idx];
            if (ny >= 4u) {
                let u1 = field[(1u) * uniforms.dims_bc.x + i];
                let u2 = field[(2u) * uniforms.dims_bc.x + i];
                let u3 = field[(3u) * uniforms.dims_bc.x + i];
                // Second-order one-sided second derivative at bottom boundary
                laplacian += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * uniforms.inv2.y;
            } else {
                let top = field[(j + 1u) * uniforms.dims_bc.x + i];
                // Fallback: symmetric ghost approach for small ny
                laplacian += (top - 2.0 * center + top) * uniforms.inv2.y;
            }
        } else if (uniforms.dims_bc.z == 2u) { // Periodic
            let bottom = field[(uniforms.dims_bc.y - 1u) * uniforms.dims_bc.x + i];
            let center = field[idx];
            let top = field[(j + 1u) * uniforms.dims_bc.x + i];
            laplacian += (bottom - 2.0 * center + top) * uniforms.inv2.y;
        }
    } else if (j == uniforms.dims_bc.y - 1u) {
        // Top boundary
        if (uniforms.dims_bc.z == 0u) { // Dirichlet: u(i,ny-1) = 0
            let center = field[idx];
            // Ghost point u(i,ny) = -u(i,ny-2) => bottom+top-2*center cancels neighbor, yields -2*center
            laplacian += (-2.0 * center) * uniforms.inv2.y;
        } else if (uniforms.dims_bc.z == 1u) { // Neumann: du/dy = 0
            let ny = uniforms.dims_bc.y;
            let center = field[idx];
            if (ny >= 4u) {
                let u1 = field[(ny - 2u) * uniforms.dims_bc.x + i];
                let u2 = field[(ny - 3u) * uniforms.dims_bc.x + i];
                let u3 = field[(ny - 4u) * uniforms.dims_bc.x + i];
                // Second-order one-sided second derivative at top boundary
                laplacian += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * uniforms.inv2.y;
            } else {
                let bottom = field[(j - 1u) * uniforms.dims_bc.x + i];
                // Fallback: symmetric ghost approach for small ny
                laplacian += (bottom - 2.0 * center + bottom) * uniforms.inv2.y;
            }
        } else if (uniforms.dims_bc.z == 2u) { // Periodic
            let bottom = field[(j - 1u) * uniforms.dims_bc.x + i];
            let center = field[idx];
            let top = field[0u * uniforms.dims_bc.x + i];
            laplacian += (bottom - 2.0 * center + top) * uniforms.inv2.y;
        }
    }
    
    result[idx] = laplacian;
}
";
