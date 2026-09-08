// 2D Laplacian operator compute shader (mathematically rigorous)
// Implements ∇²u with proper boundary conditions:
// - Dirichlet: ghost points enforce u=0 ⇒ −2u_boundary per axis
// - Neumann: second-order one-sided second derivatives (fallback for small grids)
// - Periodic: wrap neighbors across endpoints

struct Params {
    nx: u32,
    ny: u32,
    bc_x: u32,
    bc_y: u32,
    dx_inv2: f32,
    dy_inv2: f32,
    _pad0: f32,
    _pad1: f32,
};

@group(0) @binding(0)
var<uniform> params: Params;

@group(0) @binding(1)
var<storage, read> input: array<f32>;

@group(0) @binding(3)
var<storage, read_write> output: array<f32>;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;

    if (i >= params.nx || j >= params.ny) {
        return;
    }

    let idx = j * params.nx + i;
    var lap: f32 = 0.0;
    let center = input[idx];

    // X-direction contribution
    if (params.bc_x == 2u) {
        // Periodic in X — 5-point stencil across domain
        let nx = params.nx;
        let i_m1 = (i + nx - 1u) % nx;
        let i_m2 = (i + nx - 2u) % nx;
        let i_p1 = (i + 1u) % nx;
        let i_p2 = (i + 2u) % nx;
        let u_m2 = input[j * nx + i_m2];
        let u_m1 = input[j * nx + i_m1];
        let u_p1 = input[j * nx + i_p1];
        let u_p2 = input[j * nx + i_p2];
        // (-1/12*u_{i-2} + 4/3*u_{i-1} - 5/2*u_i + 4/3*u_{i+1} - 1/12*u_{i+2}) / dx^2
        let c_m2 = -1.0 / 12.0;
        let c_m1 =  4.0 /  3.0;
        let c_0  = -5.0 /  2.0;
        let c_p1 =  4.0 /  3.0;
        let c_p2 = -1.0 / 12.0;
        lap += (c_m2 * u_m2 + c_m1 * u_m1 + c_0 * center + c_p1 * u_p1 + c_p2 * u_p2) * params.dx_inv2;
    } else if (i > 1u && i < params.nx - 2u) {
        let u_m2 = input[idx - 2u];
        let u_m1 = input[idx - 1u];
        let u_p1 = input[idx + 1u];
        let u_p2 = input[idx + 2u];
        let c_m2 = -1.0 / 12.0;
        let c_m1 =  4.0 /  3.0;
        let c_0  = -5.0 /  2.0;
        let c_p1 =  4.0 /  3.0;
        let c_p2 = -1.0 / 12.0;
        lap += (c_m2 * u_m2 + c_m1 * u_m1 + c_0 * center + c_p1 * u_p1 + c_p2 * u_p2) * params.dx_inv2;
    } else if (i > 0u && i < params.nx - 1u) {
        // Interior (Dirichlet/Neumann): 3-point central difference with numerically stable form
        // Avoid catastrophic cancellation by summing first-order differences
        let left = input[idx - 1u];
        let right = input[idx + 1u];
        // Use fused multiply-add to minimize rounding in (left + right - 2*center)
        let t_lr = left + right;
        lap += fma(-2.0, center, t_lr) * params.dx_inv2;
    } else if (i == 0u) {
        if (params.bc_x == 0u) {
            // Dirichlet: outside value = 0 ⇒ (0 - 2*center + right)
            let right = input[idx + 1u];
            lap += (0.0 - 2.0 * center + right) * params.dx_inv2;
        } else { // Neumann
            // du/dx=0, second-order one-sided second derivative
            let nx = params.nx;
            if (nx >= 4u) {
                let u1 = input[j * nx + 1u];
                let u2 = input[j * nx + 2u];
                let u3 = input[j * nx + 3u];
                lap += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * params.dx_inv2;
            } else {
                let right = input[idx + 1u];
                lap += (right - 2.0 * center + right) * params.dx_inv2;
            }
        }
    } else { // i == nx-1
        if (params.bc_x == 0u) {
            // Dirichlet: outside value = 0 ⇒ (left - 2*center + 0)
            let left = input[idx - 1u];
            lap += (left - 2.0 * center + 0.0) * params.dx_inv2;
        } else { // Neumann
            let nx = params.nx;
            if (nx >= 4u) {
                let u1 = input[j * nx + (nx - 2u)];
                let u2 = input[j * nx + (nx - 3u)];
                let u3 = input[j * nx + (nx - 4u)];
                lap += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * params.dx_inv2;
            } else {
                let left = input[idx - 1u];
                lap += (left - 2.0 * center + left) * params.dx_inv2;
            }
        }
    }

    // Y-direction contribution
    if (params.bc_y == 2u) {
        // Periodic in Y — 5-point stencil across domain
        let ny = params.ny;
        let j_m1 = (j + ny - 1u) % ny;
        let j_m2 = (j + ny - 2u) % ny;
        let j_p1 = (j + 1u) % ny;
        let j_p2 = (j + 2u) % ny;
        let u_m2 = input[j_m2 * params.nx + i];
        let u_m1 = input[j_m1 * params.nx + i];
        let u_p1 = input[j_p1 * params.nx + i];
        let u_p2 = input[j_p2 * params.nx + i];
        let c_m2y = -1.0 / 12.0;
        let c_m1y =  4.0 /  3.0;
        let c_0y  = -5.0 /  2.0;
        let c_p1y =  4.0 /  3.0;
        let c_p2y = -1.0 / 12.0;
        lap += (c_m2y * u_m2 + c_m1y * u_m1 + c_0y * center + c_p1y * u_p1 + c_p2y * u_p2) * params.dy_inv2;
    } else if (j > 1u && j < params.ny - 2u) {
        let u_m2 = input[(j - 2u) * params.nx + i];
        let u_m1 = input[(j - 1u) * params.nx + i];
        let u_p1 = input[(j + 1u) * params.nx + i];
        let u_p2 = input[(j + 2u) * params.nx + i];
        let c_m2y = -1.0 / 12.0;
        let c_m1y =  4.0 /  3.0;
        let c_0y  = -5.0 /  2.0;
        let c_p1y =  4.0 /  3.0;
        let c_p2y = -1.0 / 12.0;
        lap += (c_m2y * u_m2 + c_m1y * u_m1 + c_0y * center + c_p1y * u_p1 + c_p2y * u_p2) * params.dy_inv2;
    } else if (j > 0u && j < params.ny - 1u) {
        // Interior (Dirichlet/Neumann): 3-point central difference with numerically stable form
        let bottom = input[(j - 1u) * params.nx + i];
        let top = input[(j + 1u) * params.nx + i];
        let t_bt = bottom + top;
        lap += fma(-2.0, center, t_bt) * params.dy_inv2;
    } else if (j == 0u) {
        if (params.bc_y == 0u) {
            // Dirichlet: outside value = 0 ⇒ (0 - 2*center + top)
            let top = input[(j + 1u) * params.nx + i];
            lap += (0.0 - 2.0 * center + top) * params.dy_inv2;
        } else { // Neumann
            let ny = params.ny;
            if (ny >= 4u) {
                let u1 = input[1u * params.nx + i];
                let u2 = input[2u * params.nx + i];
                let u3 = input[3u * params.nx + i];
                lap += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * params.dy_inv2;
            } else {
                let top = input[(j + 1u) * params.nx + i];
                lap += (top - 2.0 * center + top) * params.dy_inv2;
            }
        }
    } else { // j == ny-1
        if (params.bc_y == 0u) {
            // Dirichlet: outside value = 0 ⇒ (bottom - 2*center + 0)
            let bottom = input[(j - 1u) * params.nx + i];
            lap += (bottom - 2.0 * center + 0.0) * params.dy_inv2;
        } else { // Neumann
            let ny = params.ny;
            if (ny >= 4u) {
                let u1 = input[(ny - 2u) * params.nx + i];
                let u2 = input[(ny - 3u) * params.nx + i];
                let u3 = input[(ny - 4u) * params.nx + i];
                lap += (2.0 * center - 5.0 * u1 + 4.0 * u2 - u3) * params.dy_inv2;
            } else {
                let bottom = input[(j - 1u) * params.nx + i];
                lap += (bottom - 2.0 * center + bottom) * params.dy_inv2;
            }
        }
    }

    // Positive Laplacian
    output[idx] = lap;
}
