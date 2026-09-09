// Grid-cutoff branch of the detached-eddy-simulation length scale.

struct Params {
    dimensions: vec4<u32>,
    scales: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(8, 8, 1)
fn des_grid_scale(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    if (i >= params.dimensions.x || j >= params.dimensions.y) {
        return;
    }
    let index = j * params.dimensions.x + i;
    output[index] = params.scales.w * params.scales.z;
}
