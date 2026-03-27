// Pick shader — renders flat color IDs for GPU-based selection picking.
//
// Each triangle is rendered with a unique color encoding (node_index, face_index)
// as RGBA channels. The CPU reads back the pixel at the click location to
// determine which entity was selected.

struct ViewUniforms {
    view_proj: mat4x4<f32>,
    camera_pos: vec3<f32>,
    _pad0: f32,
    light_dir: vec3<f32>,
    _pad1: f32,
    shading_mode: u32,
    field_min: f32,
    field_max: f32,
    field_active: u32,
};

@group(0) @binding(0)
var<uniform> uniforms: ViewUniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) region_id: u32,
    @location(3) field_value: f32,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) @interpolate(flat) v_region_id: u32,
    @location(1) @interpolate(flat) v_face_id: u32,
};

// The face_id is encoded in the field_value as a bitcast f32 → u32.
@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.clip_position = uniforms.view_proj * vec4<f32>(input.position, 1.0);
    out.v_region_id = input.region_id;
    out.v_face_id = bitcast<u32>(input.field_value);
    return out;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    // Encode (node_index=region_id, face_index) into RGBA bytes.
    let node = input.v_region_id;
    let face = input.v_face_id;
    let r = f32(node & 0xFFu) / 255.0;
    let g = f32((node >> 8u) & 0xFFu) / 255.0;
    let b = f32(face & 0xFFu) / 255.0;
    let a = f32((face >> 8u) & 0xFFu) / 255.0;
    return vec4<f32>(r, g, b, a);
}
