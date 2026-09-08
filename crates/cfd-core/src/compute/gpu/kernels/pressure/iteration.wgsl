// Weighted-Jacobi iteration for a three-dimensional pressure Poisson equation.

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
fn pressure_iteration(@builtin(global_invocation_id) global_id: vec3<u32>) {
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
        let interior_i = clamp(i, 1u, params.dimensions.x - 2u);
        let interior_j = clamp(j, 1u, params.dimensions.y - 2u);
        let interior_k = clamp(k, 1u, params.dimensions.z - 2u);
        output[index] = pressure[index_3d(interior_i, interior_j, interior_k)];
        return;
    }

    let inverse_x = params.inverse_spacing_squared_relaxation.x;
    let inverse_y = params.inverse_spacing_squared_relaxation.y;
    let inverse_z = params.inverse_spacing_squared_relaxation.z;
    let neighbor_sum = (pressure[index_3d(i - 1u, j, k)]
        + pressure[index_3d(i + 1u, j, k)]) * inverse_x
        + (pressure[index_3d(i, j - 1u, k)]
        + pressure[index_3d(i, j + 1u, k)]) * inverse_y
        + (pressure[index_3d(i, j, k - 1u)]
        + pressure[index_3d(i, j, k + 1u)]) * inverse_z;
    let jacobi = (neighbor_sum - source[index]) / (2.0 * (inverse_x + inverse_y + inverse_z));
    let relaxation = params.inverse_spacing_squared_relaxation.w;
    output[index] = (1.0 - relaxation) * pressure[index] + relaxation * jacobi;
}
