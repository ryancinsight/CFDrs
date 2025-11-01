// 2D Momentum operator compute shader
// Computes simplified momentum equation: -ν∇²u + (u·∇)u

struct Params {
    nx: u32,
    ny: u32,
    nz: u32,  // unused for 2D
    dx: f32,
    dy: f32,
    dz: f32,  // unused for 2D
    viscosity: f32,
    density: f32,  // unused in this simplified version
};

@group(0) @binding(0)
var<uniform> params: Params;

@group(0) @binding(1)
var<storage, read> input: array<f32>;  // velocity field u

@group(0) @binding(2)
var<storage, read> velocity_u: array<f32>;  // convection velocity u

@group(0) @binding(3)
var<storage, read> velocity_v: array<f32>;  // convection velocity v

@group(0) @binding(4)
var<storage, read_write> output: array<f32>;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;

    // Check bounds
    if (i >= params.nx || j >= params.ny) {
        return;
    }

    let idx = j * params.nx + i;
    let center = input[idx];

    // Compute derivatives
    let dx2_inv = 1.0 / (params.dx * params.dx);
    let dy2_inv = 1.0 / (params.dy * params.dy);
    let dx_inv = 1.0 / params.dx;
    let dy_inv = 1.0 / params.dy;

    var laplacian: f32 = 0.0;
    var convection: f32 = 0.0;

    // Laplacian term: ∇²u
    if (i > 0u && i < params.nx - 1u && j > 0u && j < params.ny - 1u) {
        // Interior points
        let d2u_dx2 = (input[idx - 1u] + input[idx + 1u] - 2.0 * center) * dx2_inv;
        let d2u_dy2 = (input[idx - params.nx] + input[idx + params.nx] - 2.0 * center) * dy2_inv;
        laplacian = d2u_dx2 + d2u_dy2;

        // Convection term: (u·∇)u
        let ui = velocity_u[idx];
        let vi = velocity_v[idx];

        // Central differences for derivatives
        let du_dx = (input[idx + 1u] - input[idx - 1u]) * 0.5 * dx_inv;
        let du_dy = (input[idx + params.nx] - input[idx - params.nx]) * 0.5 * dy_inv;

        convection = ui * du_dx + vi * du_dy;
    } else {
        // Boundary conditions (no-slip)
        output[idx] = 0.0;
        return;
    }

    // Combine terms: result = -ν∇²u + convection
    output[idx] = -params.viscosity * laplacian + convection;
}
