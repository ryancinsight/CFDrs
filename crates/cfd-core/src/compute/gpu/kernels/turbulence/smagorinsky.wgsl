// Two-dimensional Smagorinsky subgrid viscosity.

struct Params {
    dimensions: vec4<u32>,
    scales: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> velocity_x: array<f32>;
@group(0) @binding(2) var<storage, read> velocity_y: array<f32>;
@group(0) @binding(3) var<storage, read_write> output: array<f32>;

fn index_2d(i: u32, j: u32) -> u32 {
    return j * params.dimensions.x + i;
}

@compute @workgroup_size(8, 8, 1)
fn smagorinsky_viscosity(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    if (i >= params.dimensions.x || j >= params.dimensions.y) {
        return;
    }
    let index = index_2d(i, j);
    if (i == 0u || i == params.dimensions.x - 1u
        || j == 0u || j == params.dimensions.y - 1u) {
        output[index] = 0.0;
        return;
    }

    let derivative_x_x = (velocity_x[index_2d(i + 1u, j)]
        - velocity_x[index_2d(i - 1u, j)]) * params.scales.x;
    let derivative_x_y = (velocity_x[index_2d(i, j + 1u)]
        - velocity_x[index_2d(i, j - 1u)]) * params.scales.y;
    let derivative_y_x = (velocity_y[index_2d(i + 1u, j)]
        - velocity_y[index_2d(i - 1u, j)]) * params.scales.x;
    let derivative_y_y = (velocity_y[index_2d(i, j + 1u)]
        - velocity_y[index_2d(i, j - 1u)]) * params.scales.y;
    let shear = 0.5 * (derivative_x_y + derivative_y_x);
    let strain_magnitude = sqrt(2.0 * (derivative_x_x * derivative_x_x
        + derivative_y_y * derivative_y_y + 2.0 * shear * shear));
    let filter_factor = params.scales.w * params.scales.z;
    output[index] = filter_factor * filter_factor * strain_magnitude;
}
