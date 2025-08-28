// Diffusion kernel for CFD - Central difference scheme
// Solves: ∂u/∂t = ν∇²u

struct DiffusionParams {
    nx: u32,
    ny: u32,
    nz: u32,
    dx: f32,
    dy: f32,
    dz: f32,
    dt: f32,
    nu: f32,  // Kinematic viscosity
}

@group(0) @binding(0) var<uniform> params: DiffusionParams;
@group(0) @binding(1) var<storage, read> u_in: array<f32>;
@group(0) @binding(2) var<storage, read_write> u_out: array<f32>;

fn idx_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.nx * params.ny + j * params.nx + i;
}

@compute @workgroup_size(8, 8, 1)
fn diffusion_central(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    
    // Boundary check
    if (i >= params.nx || j >= params.ny || k >= params.nz) {
        return;
    }
    
    // Skip boundaries (Dirichlet BC)
    if (i == 0u || i == params.nx - 1u || 
        j == 0u || j == params.ny - 1u ||
        k == 0u || k == params.nz - 1u) {
        let idx = idx_3d(i, j, k);
        u_out[idx] = u_in[idx];
        return;
    }
    
    let idx = idx_3d(i, j, k);
    let idx_ip = idx_3d(i + 1u, j, k);
    let idx_im = idx_3d(i - 1u, j, k);
    let idx_jp = idx_3d(i, j + 1u, k);
    let idx_jm = idx_3d(i, j - 1u, k);
    let idx_kp = idx_3d(i, j, k + 1u);
    let idx_km = idx_3d(i, j, k - 1u);
    
    // Central difference for Laplacian
    let d2u_dx2 = (u_in[idx_ip] - 2.0 * u_in[idx] + u_in[idx_im]) / (params.dx * params.dx);
    let d2u_dy2 = (u_in[idx_jp] - 2.0 * u_in[idx] + u_in[idx_jm]) / (params.dy * params.dy);
    let d2u_dz2 = (u_in[idx_kp] - 2.0 * u_in[idx] + u_in[idx_km]) / (params.dz * params.dz);
    
    let laplacian = d2u_dx2 + d2u_dy2 + d2u_dz2;
    
    // Explicit time stepping
    u_out[idx] = u_in[idx] + params.dt * params.nu * laplacian;
}