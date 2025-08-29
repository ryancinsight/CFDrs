// Velocity correction kernel for CFD - SIMPLE algorithm
// Updates velocity based on pressure gradient
// Reference: Patankar (1980) - Numerical Heat Transfer and Fluid Flow

struct GridParams {
    nx: u32,
    ny: u32,
    nz: u32,
    dx: f32,
    dy: f32,
    dz: f32,
    dt: f32,
    density: f32,
}

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<storage, read> u_star: array<f32>;  // Intermediate velocity x
@group(0) @binding(2) var<storage, read> v_star: array<f32>;  // Intermediate velocity y
@group(0) @binding(3) var<storage, read> w_star: array<f32>;  // Intermediate velocity z
@group(0) @binding(4) var<storage, read> pressure: array<f32>;
@group(0) @binding(5) var<storage, read_write> u_out: array<f32>;
@group(0) @binding(6) var<storage, read_write> v_out: array<f32>;
@group(0) @binding(7) var<storage, read_write> w_out: array<f32>;

fn idx_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.nx * params.ny + j * params.nx + i;
}

@compute @workgroup_size(8, 8, 1)
fn velocity_correction(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    
    // Boundary check
    if (i >= params.nx || j >= params.ny || k >= params.nz) {
        return;
    }
    
    let idx = idx_3d(i, j, k);
    
    // Apply boundary conditions
    if (i == 0u || i == params.nx - 1u || 
        j == 0u || j == params.ny - 1u ||
        k == 0u || k == params.nz - 1u) {
        // Wall boundary conditions (no-slip)
        u_out[idx] = 0.0;
        v_out[idx] = 0.0;
        w_out[idx] = 0.0;
        return;
    }
    
    // Interior points - apply pressure correction
    // u^(n+1) = u* - dt/ρ * ∂p/∂x
    // v^(n+1) = v* - dt/ρ * ∂p/∂y
    // w^(n+1) = w* - dt/ρ * ∂p/∂z
    
    // Pressure gradients (central difference)
    let idx_im = idx_3d(i - 1u, j, k);
    let idx_ip = idx_3d(i + 1u, j, k);
    let dp_dx = (pressure[idx_ip] - pressure[idx_im]) / (2.0 * params.dx);
    
    let idx_jm = idx_3d(i, j - 1u, k);
    let idx_jp = idx_3d(i, j + 1u, k);
    let dp_dy = (pressure[idx_jp] - pressure[idx_jm]) / (2.0 * params.dy);
    
    let idx_km = idx_3d(i, j, k - 1u);
    let idx_kp = idx_3d(i, j, k + 1u);
    let dp_dz = (pressure[idx_kp] - pressure[idx_km]) / (2.0 * params.dz);
    
    // Apply velocity correction
    let correction_factor = params.dt / params.density;
    u_out[idx] = u_star[idx] - correction_factor * dp_dx;
    v_out[idx] = v_star[idx] - correction_factor * dp_dy;
    w_out[idx] = w_star[idx] - correction_factor * dp_dz;
}

@compute @workgroup_size(8, 8, 1)
fn divergence_compute(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    
    // Boundary check
    if (i >= params.nx || j >= params.ny || k >= params.nz) {
        return;
    }
    
    let idx = idx_3d(i, j, k);
    
    // Boundaries have zero divergence
    if (i == 0u || i == params.nx - 1u || 
        j == 0u || j == params.ny - 1u ||
        k == 0u || k == params.nz - 1u) {
        u_out[idx] = 0.0;  // Using u_out to store divergence
        return;
    }
    
    // Compute divergence: ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
    let idx_im = idx_3d(i - 1u, j, k);
    let idx_ip = idx_3d(i + 1u, j, k);
    let du_dx = (u_star[idx_ip] - u_star[idx_im]) / (2.0 * params.dx);
    
    let idx_jm = idx_3d(i, j - 1u, k);
    let idx_jp = idx_3d(i, j + 1u, k);
    let dv_dy = (v_star[idx_jp] - v_star[idx_jm]) / (2.0 * params.dy);
    
    let idx_km = idx_3d(i, j, k - 1u);
    let idx_kp = idx_3d(i, j, k + 1u);
    let dw_dz = (w_star[idx_kp] - w_star[idx_km]) / (2.0 * params.dz);
    
    // Store divergence (will be used as source term for pressure Poisson equation)
    u_out[idx] = params.density * (du_dx + dv_dy + dw_dz) / params.dt;
}