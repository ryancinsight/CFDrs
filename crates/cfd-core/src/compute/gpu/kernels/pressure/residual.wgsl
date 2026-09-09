// Absolute pointwise residual for a three-dimensional pressure Poisson equation.

struct Params {
    dimensions: vec4<u32>,
    inverse_spacing_squared_relaxation: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> pressure: array<f32>;
@group(0) @binding(2) var<storage, read> source: array<f32>;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

fn index_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.dimensions.x * params.dimensions.y + j * params.dimensions.x + i;
}

@compute @workgroup_size(8, 8, 1)
fn pressure_residual(@builtin(global_invocation_id) global_id: vec3<u32>) {
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
        output[index] = 0.0;
        return;
    }

    let center = pressure[index];
    let second_x = (pressure[index_3d(i - 1u, j, k)] - 2.0 * center
        + pressure[index_3d(i + 1u, j, k)]) * params.inverse_spacing_squared_relaxation.x;
    let second_y = (pressure[index_3d(i, j - 1u, k)] - 2.0 * center
        + pressure[index_3d(i, j + 1u, k)]) * params.inverse_spacing_squared_relaxation.y;
    let second_z = (pressure[index_3d(i, j, k - 1u)] - 2.0 * center
        + pressure[index_3d(i, j, k + 1u)]) * params.inverse_spacing_squared_relaxation.z;
    output[index] = abs(second_x + second_y + second_z - source[index]);
}
