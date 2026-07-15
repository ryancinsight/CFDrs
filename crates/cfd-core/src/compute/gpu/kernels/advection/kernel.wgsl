// First-order upwind advection on independent z-planes.

struct Params {
    dimensions: vec4<u32>,
    spacing_step: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> scalar: array<f32>;
@group(0) @binding(2) var<storage, read> velocity_x: array<f32>;
@group(0) @binding(3) var<storage, read> velocity_y: array<f32>;
@group(0) @binding(4) var<storage, read_write> output: array<f32>;

fn index_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.dimensions.x * params.dimensions.y + j * params.dimensions.x + i;
}

@compute @workgroup_size(8, 8, 1)
fn advection_upwind(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    if (i >= params.dimensions.x || j >= params.dimensions.y || k >= params.dimensions.z) {
        return;
    }

    let index = index_3d(i, j, k);
    if (i == 0u || i == params.dimensions.x - 1u || j == 0u || j == params.dimensions.y - 1u) {
        output[index] = scalar[index];
        return;
    }

    let vx = velocity_x[index];
    let vy = velocity_y[index];
    var derivative_x: f32;
    var derivative_y: f32;
    if (vx > 0.0) {
        derivative_x = (scalar[index] - scalar[index_3d(i - 1u, j, k)]) / params.spacing_step.x;
    } else {
        derivative_x = (scalar[index_3d(i + 1u, j, k)] - scalar[index]) / params.spacing_step.x;
    }
    if (vy > 0.0) {
        derivative_y = (scalar[index] - scalar[index_3d(i, j - 1u, k)]) / params.spacing_step.y;
    } else {
        derivative_y = (scalar[index_3d(i, j + 1u, k)] - scalar[index]) / params.spacing_step.y;
    }

    output[index] = scalar[index]
        - params.spacing_step.w * (vx * derivative_x + vy * derivative_y);
}
