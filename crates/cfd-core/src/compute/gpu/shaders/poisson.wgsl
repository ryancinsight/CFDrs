// Poisson equation solver using Jacobi iteration
// Solves ∇²φ = f

struct PoissonParams {
    nx: u32,
    ny: u32,
    dx: f32,
    dy: f32,
    omega: f32, // SOR relaxation parameter
}

@group(0) @binding(0)
var<uniform> params: PoissonParams;

@group(0) @binding(1)
var<storage, read> phi_in: array<f32>;

@group(0) @binding(2)
var<storage, read> source: array<f32>;

@group(0) @binding(3)
var<storage, read_write> phi_out: array<f32>;

fn get_index(i: u32, j: u32) -> u32 {
    return i * params.ny + j;
}

@compute @workgroup_size(8, 8)
fn jacobi_iteration(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    
    // Boundary check
    if (i == 0u || i >= params.nx - 1u || j == 0u || j >= params.ny - 1u) {
        // Copy boundary values unchanged
        if (i < params.nx && j < params.ny) {
            let idx = get_index(i, j);
            phi_out[idx] = phi_in[idx];
        }
        return;
    }
    
    // Interior point calculation
    let idx = get_index(i, j);
    let idx_left = get_index(i - 1u, j);
    let idx_right = get_index(i + 1u, j);
    let idx_bottom = get_index(i, j - 1u);
    let idx_top = get_index(i, j + 1u);
    
    let dx2 = params.dx * params.dx;
    let dy2 = params.dy * params.dy;
    let factor = 0.5 / (1.0 / dx2 + 1.0 / dy2);
    
    // Jacobi update
    let new_value = factor * (
        (phi_in[idx_left] + phi_in[idx_right]) / dx2 +
        (phi_in[idx_bottom] + phi_in[idx_top]) / dy2 -
        source[idx]
    );
    
    // SOR relaxation
    phi_out[idx] = (1.0 - params.omega) * phi_in[idx] + params.omega * new_value;
}

// Red-Black Gauss-Seidel for better convergence
@compute @workgroup_size(8, 8)
fn red_black_iteration(
    @builtin(global_invocation_id) global_id: vec3<u32>,
    @builtin(workgroup_id) workgroup_id: vec3<u32>
) {
    let i = global_id.x;
    let j = global_id.y;
    
    // Determine if this is a red or black point
    let is_red = (i + j) % 2u == 0u;
    
    // Skip if wrong color for this pass
    // (controlled by dispatch parameter)
    
    if (i == 0u || i >= params.nx - 1u || j == 0u || j >= params.ny - 1u) {
        return;
    }
    
    let idx = get_index(i, j);
    let idx_left = get_index(i - 1u, j);
    let idx_right = get_index(i + 1u, j);
    let idx_bottom = get_index(i, j - 1u);
    let idx_top = get_index(i, j + 1u);
    
    let dx2 = params.dx * params.dx;
    let dy2 = params.dy * params.dy;
    let factor = params.omega / (2.0 * (1.0 / dx2 + 1.0 / dy2));
    
    let residual = (phi_in[idx_left] + phi_in[idx_right]) / dx2 +
                   (phi_in[idx_bottom] + phi_in[idx_top]) / dy2 -
                   source[idx];
    
    phi_out[idx] = (1.0 - params.omega) * phi_in[idx] + factor * residual;
}

// Calculate residual for convergence check
@compute @workgroup_size(8, 8)
fn calculate_residual(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    
    if (i == 0u || i >= params.nx - 1u || j == 0u || j >= params.ny - 1u) {
        return;
    }
    
    let idx = get_index(i, j);
    let idx_left = get_index(i - 1u, j);
    let idx_right = get_index(i + 1u, j);
    let idx_bottom = get_index(i, j - 1u);
    let idx_top = get_index(i, j + 1u);
    
    let dx2 = params.dx * params.dx;
    let dy2 = params.dy * params.dy;
    
    // Calculate Laplacian
    let laplacian = (phi_in[idx_left] - 2.0 * phi_in[idx] + phi_in[idx_right]) / dx2 +
                    (phi_in[idx_bottom] - 2.0 * phi_in[idx] + phi_in[idx_top]) / dy2;
    
    // Residual = ∇²φ - f
    phi_out[idx] = abs(laplacian - source[idx]);
}