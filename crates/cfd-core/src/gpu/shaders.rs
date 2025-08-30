//! GPU shader definitions for field operations

/// WGSL shader for field addition
pub const FIELD_ADD_SHADER: &str = r#"
@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn add_fields(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&a)) {
        return;
    }
    result[index] = a[index] + b[index];
}
"#;

/// WGSL shader for scalar multiplication
pub const FIELD_MUL_SHADER: &str = r#"
@group(0) @binding(0) var<storage, read> field: array<f32>;
@group(0) @binding(1) var<uniform> scalar: f32;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn multiply_field(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&field)) {
        return;
    }
    result[index] = field[index] * scalar;
}
"#;

/// WGSL shader for 2D Laplacian
pub const LAPLACIAN_2D_SHADER: &str = r#"
struct Uniforms {
    nx: u32,
    ny: u32,
    dx_inv2: f32,
    dy_inv2: f32,
}

@group(0) @binding(0) var<storage, read> field: array<f32>;
@group(0) @binding(1) var<uniform> uniforms: Uniforms;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(8, 8)
fn laplacian_2d(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    
    if (i >= uniforms.nx || j >= uniforms.ny) {
        return;
    }
    
    let idx = j * uniforms.nx + i;
    var laplacian = 0.0;
    
    // X direction
    if (i > 0u && i < uniforms.nx - 1u) {
        let left = field[j * uniforms.nx + (i - 1u)];
        let center = field[idx];
        let right = field[j * uniforms.nx + (i + 1u)];
        laplacian += (left - 2.0 * center + right) * uniforms.dx_inv2;
    }
    
    // Y direction
    if (j > 0u && j < uniforms.ny - 1u) {
        let bottom = field[(j - 1u) * uniforms.nx + i];
        let center = field[idx];
        let top = field[(j + 1u) * uniforms.nx + i];
        laplacian += (bottom - 2.0 * center + top) * uniforms.dy_inv2;
    }
    
    result[idx] = laplacian;
}
"#;