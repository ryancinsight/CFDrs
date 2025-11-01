// 2D Laplacian operator compute shader
// Computes -∇²x where ∇² is the discrete Laplacian operator

struct Params {
    nx: u32,
    ny: u32,
    nz: u32,  // unused for 2D
    dx: f32,
    dy: f32,
    dz: f32,  // unused for 2D
    dt: f32,  // unused
};

@group(0) @binding(0)
var<uniform> params: Params;

@group(0) @binding(1)
var<storage, read> input: array<f32>;

@group(0) @binding(3)
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

    // Compute second derivatives
    let dx2_inv = 1.0 / (params.dx * params.dx);
    let dy2_inv = 1.0 / (params.dy * params.dy);

    var d2x_dx2: f32 = 0.0;
    var d2x_dy2: f32 = 0.0;

    // Second derivative in x-direction
    if (i > 0u && i < params.nx - 1u) {
        d2x_dx2 = (input[idx - 1u] + input[idx + 1u] - 2.0 * input[idx]) * dx2_inv;
    } else {
        // Neumann boundary conditions
        if (i == 0u && params.nx > 1u) {
            d2x_dx2 = (input[idx + 1u] + input[idx + 1u] - 2.0 * input[idx]) * dx2_inv;
        } else if (i == params.nx - 1u && params.nx > 1u) {
            d2x_dx2 = (input[idx - 1u] + input[idx - 1u] - 2.0 * input[idx]) * dx2_inv;
        }
    }

    // Second derivative in y-direction
    if (j > 0u && j < params.ny - 1u) {
        let idx_north = (j + 1u) * params.nx + i;
        let idx_south = (j - 1u) * params.nx + i;
        d2x_dy2 = (input[idx_north] + input[idx_south] - 2.0 * input[idx]) * dy2_inv;
    } else {
        // Neumann boundary conditions
        if (j == 0u && params.ny > 1u) {
            let idx_north = (j + 1u) * params.nx + i;
            d2x_dy2 = (input[idx_north] + input[idx_north] - 2.0 * input[idx]) * dy2_inv;
        } else if (j == params.ny - 1u && params.ny > 1u) {
            let idx_south = (j - 1u) * params.nx + i;
            d2x_dy2 = (input[idx_south] + input[idx_south] - 2.0 * input[idx]) * dy2_inv;
        }
    }

    // Negative Laplacian: -∇²x = - (d²x/dx² + d²x/dy²)
    output[idx] = -(d2x_dx2 + d2x_dy2);
}
