// Advection kernel for CFD - Upwind scheme
// Solves: ∂u/∂t + v·∇u = 0

struct GridParams {
    nx: u32,
    ny: u32,
    nz: u32,
    dx: f32,
    dy: f32,
    dz: f32,
    dt: f32,
}

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<storage, read> u_in: array<f32>;
@group(0) @binding(2) var<storage, read> velocity_x: array<f32>;
@group(0) @binding(3) var<storage, read> velocity_y: array<f32>;
@group(0) @binding(4) var<storage, read_write> u_out: array<f32>;

fn idx_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.nx * params.ny + j * params.nx + i;
}

@compute @workgroup_size(8, 8, 1)
fn advection_upwind(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    
    // Boundary check
    if (i >= params.nx || j >= params.ny || k >= params.nz) {
        return;
    }
    
    // Skip boundaries
    if (i == 0u || i == params.nx - 1u || j == 0u || j == params.ny - 1u) {
        let idx = idx_3d(i, j, k);
        u_out[idx] = u_in[idx];
        return;
    }
    
    let idx = idx_3d(i, j, k);
    let vx = velocity_x[idx];
    let vy = velocity_y[idx];
    
    var du_dx: f32;
    var du_dy: f32;
    
    // Upwind scheme based on velocity direction
    if (vx > 0.0) {
        // Backward difference
        let idx_im = idx_3d(i - 1u, j, k);
        du_dx = (u_in[idx] - u_in[idx_im]) / params.dx;
    } else {
        // Forward difference
        let idx_ip = idx_3d(i + 1u, j, k);
        du_dx = (u_in[idx_ip] - u_in[idx]) / params.dx;
    }
    
    if (vy > 0.0) {
        // Backward difference
        let idx_jm = idx_3d(i, j - 1u, k);
        du_dy = (u_in[idx] - u_in[idx_jm]) / params.dy;
    } else {
        // Forward difference
        let idx_jp = idx_3d(i, j + 1u, k);
        du_dy = (u_in[idx_jp] - u_in[idx]) / params.dy;
    }
    
    // Update using explicit time stepping
    u_out[idx] = u_in[idx] - params.dt * (vx * du_dx + vy * du_dy);
}