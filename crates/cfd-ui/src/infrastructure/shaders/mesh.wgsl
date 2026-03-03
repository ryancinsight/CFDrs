// Mesh vertex/fragment shader — Phong illumination with per-vertex normals.
//
// Supports three shading modes via uniform:
//   0 = flat (face normal from cross product)
//   1 = smooth (interpolated vertex normal)
//   2 = wireframe (solid color, no lighting)

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
    @location(0) world_position: vec3<f32>,
    @location(1) world_normal: vec3<f32>,
    @location(2) v_region_id: u32,
    @location(3) v_field_value: f32,
};

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let world_pos = vec4<f32>(input.position, 1.0);
    out.clip_position = uniforms.view_proj * world_pos;
    out.world_position = input.position;
    out.world_normal = input.normal;
    out.v_region_id = input.region_id;
    out.v_field_value = input.field_value;
    return out;
}

// Qualitative region palette (10 colors matching the Rust-side palette).
fn region_color(id: u32) -> vec3<f32> {
    let idx = id % 10u;
    switch (idx) {
        case 0u: { return vec3<f32>(0.122, 0.467, 0.706); }
        case 1u: { return vec3<f32>(1.0, 0.498, 0.055); }
        case 2u: { return vec3<f32>(0.173, 0.627, 0.173); }
        case 3u: { return vec3<f32>(0.839, 0.153, 0.157); }
        case 4u: { return vec3<f32>(0.580, 0.404, 0.741); }
        case 5u: { return vec3<f32>(0.549, 0.337, 0.294); }
        case 6u: { return vec3<f32>(0.890, 0.467, 0.761); }
        case 7u: { return vec3<f32>(0.498, 0.498, 0.498); }
        case 8u: { return vec3<f32>(0.737, 0.741, 0.133); }
        case 9u: { return vec3<f32>(0.090, 0.745, 0.812); }
        default: { return vec3<f32>(0.5, 0.5, 0.5); }
    }
}

// Blue-Red diverging colormap for field visualization.
fn field_colormap(t: f32) -> vec3<f32> {
    let r = t;
    let b = 1.0 - t;
    let g = 1.0 - abs(2.0 * t - 1.0);
    return vec3<f32>(r, g, b);
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    // Wireframe mode: solid gray, no lighting.
    if uniforms.shading_mode == 2u {
        return vec4<f32>(0.3, 0.3, 0.3, 1.0);
    }

    // Determine base color from region or field.
    var base_color: vec3<f32>;
    if uniforms.field_active == 1u {
        let range = uniforms.field_max - uniforms.field_min;
        var t = 0.5;
        if range > 1e-8 {
            t = clamp((input.v_field_value - uniforms.field_min) / range, 0.0, 1.0);
        }
        base_color = field_colormap(t);
    } else {
        base_color = region_color(input.v_region_id);
    }

    // Compute surface normal.
    var normal: vec3<f32>;
    if uniforms.shading_mode == 0u {
        // Flat shading: compute face normal from derivatives.
        let dpdx = dpdx(input.world_position);
        let dpdy = dpdy(input.world_position);
        normal = normalize(cross(dpdx, dpdy));
    } else {
        // Smooth shading: use interpolated vertex normal.
        normal = normalize(input.world_normal);
    }

    // Simple Phong illumination.
    let ambient = 0.15;
    let light = normalize(uniforms.light_dir);
    let diffuse = max(dot(normal, light), 0.0) * 0.65;

    let view_dir = normalize(uniforms.camera_pos - input.world_position);
    let reflect_dir = reflect(-light, normal);
    let specular = pow(max(dot(view_dir, reflect_dir), 0.0), 32.0) * 0.2;

    let color = base_color * (ambient + diffuse) + vec3<f32>(specular);
    return vec4<f32>(clamp(color, vec3<f32>(0.0), vec3<f32>(1.0)), 1.0);
}
