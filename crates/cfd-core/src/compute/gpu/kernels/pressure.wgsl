// Pressure Poisson solver kernel for CFD - Jacobi iteration
// Solves: ∇²p = f (where f is the divergence of velocity)
// Reference: Patankar (1980) - SIMPLE algorithm

struct GridParams {
    nx: u32,
    ny: u32,
    nz: u32,
    dx: f32,
    dy: f32,
    dz: f32,
    omega: f32,  // Relaxation factor (typically 1.0-1.8)
}

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<storage, read> p_in: array<f32>;
@group(0) @binding(2) var<storage, read> divergence: array<f32>;
@group(0) @binding(3) var<storage, read_write> p_out: array<f32>;

fn idx_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.nx * params.ny + j * params.nx + i;
}

@compute @workgroup_size(8, 8, 1)
fn pressure_jacobi(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    
    // Boundary check
    if (i >= params.nx || j >= params.ny || k >= params.nz) {
        return;
    }
    
    // Apply Neumann boundary conditions at domain boundaries
    if (i == 0u || i == params.nx - 1u || 
        j == 0u || j == params.ny - 1u ||
        k == 0u || k == params.nz - 1u) {
        let idx = idx_3d(i, j, k);
        
        // Copy pressure from nearest interior cell (zero gradient)
        if (i == 0u) {
            p_out[idx] = p_in[idx_3d(1u, j, k)];
        } else if (i == params.nx - 1u) {
            p_out[idx] = p_in[idx_3d(params.nx - 2u, j, k)];
        } else if (j == 0u) {
            p_out[idx] = p_in[idx_3d(i, 1u, k)];
        } else if (j == params.ny - 1u) {
            p_out[idx] = p_in[idx_3d(i, params.ny - 2u, k)];
        } else if (k == 0u) {
            p_out[idx] = p_in[idx_3d(i, j, 1u)];
        } else if (k == params.nz - 1u) {
            p_out[idx] = p_in[idx_3d(i, j, params.nz - 2u)];
        } else {
            p_out[idx] = p_in[idx];
        }
        return;
    }
    
    // Interior points - Jacobi iteration
    let idx = idx_3d(i, j, k);
    let idx_im = idx_3d(i - 1u, j, k);
    let idx_ip = idx_3d(i + 1u, j, k);
    let idx_jm = idx_3d(i, j - 1u, k);
    let idx_jp = idx_3d(i, j + 1u, k);
    let idx_km = idx_3d(i, j, k - 1u);
    let idx_kp = idx_3d(i, j, k + 1u);
    
    let dx2 = params.dx * params.dx;
    let dy2 = params.dy * params.dy;
    let dz2 = params.dz * params.dz;
    
    // 3D Poisson equation discretization
    let p_jacobi = (
        (p_in[idx_im] + p_in[idx_ip]) / dx2 +
        (p_in[idx_jm] + p_in[idx_jp]) / dy2 +
        (p_in[idx_km] + p_in[idx_kp]) / dz2 -
        divergence[idx]
    ) / (2.0 / dx2 + 2.0 / dy2 + 2.0 / dz2);
    
    // Apply relaxation (SOR - Successive Over-Relaxation)
    p_out[idx] = (1.0 - params.omega) * p_in[idx] + params.omega * p_jacobi;
}

@compute @workgroup_size(256, 1, 1)
fn pressure_residual(@builtin(global_invocation_id) global_id: vec3<u32>) {
    // Compute residual for convergence check
    let idx = global_id.x;
    
    if (idx >= params.nx * params.ny * params.nz) {
        return;
    }
    
    // Convert linear index to 3D
    let k = idx / (params.nx * params.ny);
    let j = (idx % (params.nx * params.ny)) / params.nx;
    let i = idx % params.nx;
    
    // Skip boundaries
    if (i == 0u || i == params.nx - 1u || 
        j == 0u || j == params.ny - 1u ||
        k == 0u || k == params.nz - 1u) {
        p_out[idx] = 0.0;
        return;
    }
    
    let idx_im = idx_3d(i - 1u, j, k);
    let idx_ip = idx_3d(i + 1u, j, k);
    let idx_jm = idx_3d(i, j - 1u, k);
    let idx_jp = idx_3d(i, j + 1u, k);
    let idx_km = idx_3d(i, j, k - 1u);
    let idx_kp = idx_3d(i, j, k + 1u);
    
    let dx2 = params.dx * params.dx;
    let dy2 = params.dy * params.dy;
    let dz2 = params.dz * params.dz;
    
    // Compute Laplacian
    let laplacian = 
        (p_in[idx_im] - 2.0 * p_in[idx] + p_in[idx_ip]) / dx2 +
        (p_in[idx_jm] - 2.0 * p_in[idx] + p_in[idx_jp]) / dy2 +
        (p_in[idx_km] - 2.0 * p_in[idx] + p_in[idx_kp]) / dz2;
    
    // Residual = ∇²p - f
    p_out[idx] = abs(laplacian - divergence[idx]);
}