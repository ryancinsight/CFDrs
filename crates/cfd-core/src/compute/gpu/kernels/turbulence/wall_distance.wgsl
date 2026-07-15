// Exact distance to the nearest wall of a rectangular two-dimensional domain.

struct Params {
    dimensions: vec4<u32>,
    scales: vec4<f32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read_write> output: array<f32>;

@compute @workgroup_size(8, 8, 1)
fn wall_distance(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    if (i >= params.dimensions.x || j >= params.dimensions.y) {
        return;
    }
    let spacing_x = 0.5 / params.scales.x;
    let spacing_y = 0.5 / params.scales.y;
    let distance_x = min(f32(i), f32(params.dimensions.x - 1u - i)) * spacing_x;
    let distance_y = min(f32(j), f32(params.dimensions.y - 1u - j)) * spacing_y;
    output[j * params.dimensions.x + i] = min(distance_x, distance_y);
}
