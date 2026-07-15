// Centered SIMPLE velocity correction on a three-dimensional grid.

struct Params {
    dimensions: vec4<u32>,
    gradient_scale: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> velocity_x: array<f32>;
@group(0) @binding(2) var<storage, read> velocity_y: array<f32>;
@group(0) @binding(3) var<storage, read> velocity_z: array<f32>;
@group(0) @binding(4) var<storage, read> pressure: array<f32>;
@group(0) @binding(5) var<storage, read_write> output_x: array<f32>;
@group(0) @binding(6) var<storage, read_write> output_y: array<f32>;
@group(0) @binding(7) var<storage, read_write> output_z: array<f32>;

fn index_3d(i: u32, j: u32, k: u32) -> u32 {
    return k * params.dimensions.x * params.dimensions.y + j * params.dimensions.x + i;
}

@compute @workgroup_size(8, 8, 1)
fn velocity_correction(@builtin(global_invocation_id) global_id: vec3<u32>) {
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
        output_x[index] = 0.0;
        output_y[index] = 0.0;
        output_z[index] = 0.0;
        return;
    }

    let pressure_gradient_x = (pressure[index_3d(i + 1u, j, k)]
        - pressure[index_3d(i - 1u, j, k)]) * params.gradient_scale.x;
    let pressure_gradient_y = (pressure[index_3d(i, j + 1u, k)]
        - pressure[index_3d(i, j - 1u, k)]) * params.gradient_scale.y;
    let pressure_gradient_z = (pressure[index_3d(i, j, k + 1u)]
        - pressure[index_3d(i, j, k - 1u)]) * params.gradient_scale.z;
    output_x[index] = velocity_x[index] - params.gradient_scale.w * pressure_gradient_x;
    output_y[index] = velocity_y[index] - params.gradient_scale.w * pressure_gradient_y;
    output_z[index] = velocity_z[index] - params.gradient_scale.w * pressure_gradient_z;
}
