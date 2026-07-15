// Forward-Euler central-difference diffusion on a three-dimensional grid.

struct Params {
    dimensions: vec4<u32>,
    inverse_spacing_squared_step: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> input: array<f32>;
@group(0) @binding(2) var<storage, read_write> output: array<f32>;

fn index_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.dimensions.x * params.dimensions.y + j * params.dimensions.x + i;
}

@compute @workgroup_size(8, 8, 1)
fn diffusion_central(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let k = global_id.z;
    if (i >= params.dimensions.x || j >= params.dimensions.y || k >= params.dimensions.z) {
        return;
    }

    let index = index_3d(i, j, k);
    if (i == 0u || i == params.dimensions.x - 1u
        || j == 0u || j == params.dimensions.y - 1u
        || k == 0u || k == params.dimensions.z - 1u) {
        output[index] = input[index];
        return;
    }

    let center = input[index];
    let second_x = input[index_3d(i + 1u, j, k)] - 2.0 * center
        + input[index_3d(i - 1u, j, k)];
    let second_y = input[index_3d(i, j + 1u, k)] - 2.0 * center
        + input[index_3d(i, j - 1u, k)];
    let second_z = input[index_3d(i, j, k + 1u)] - 2.0 * center
        + input[index_3d(i, j, k - 1u)];
    let laplacian = second_x * params.inverse_spacing_squared_step.x
        + second_y * params.inverse_spacing_squared_step.y
        + second_z * params.inverse_spacing_squared_step.z;
    output[index] = center + params.inverse_spacing_squared_step.w * laplacian;
}
