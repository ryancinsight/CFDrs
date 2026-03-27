// Sketch vertex/fragment shader — 2D vertex-colored lines on a work plane.
//
// Transforms 2D sketch coordinates through model (work plane → world)
// and view_proj (world → clip). Supports solid and dashed line rendering.

struct SketchUniforms {
    view_proj: mat4x4<f32>,
    model: mat4x4<f32>,
    viewport_size: vec2<f32>,
    line_width_px: f32,
    _pad: f32,
};

@group(0) @binding(0)
var<uniform> uniforms: SketchUniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) color: vec4<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) v_color: vec4<f32>,
    @location(1) v_world_position: vec3<f32>,
};

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let world_pos = uniforms.model * vec4<f32>(input.position, 1.0);
    out.clip_position = uniforms.view_proj * world_pos;
    // Slight depth bias toward camera so sketch lines appear above the face.
    out.clip_position.z = out.clip_position.z - 0.0001;
    out.v_color = input.color;
    out.v_world_position = world_pos.xyz;
    return out;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    return input.v_color;
}

// Dashed variant: discard fragments in gaps.
@fragment
fn fs_dashed(input: VertexOutput) -> @location(0) vec4<f32> {
    let dash_period = 0.5;
    let dash_ratio = 0.5;
    let along = fract(length(input.v_world_position) / dash_period);
    if (along > dash_ratio) {
        discard;
    }
    return input.v_color;
}
