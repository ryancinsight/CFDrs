use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

#[derive(Clone, Copy, Debug)]
pub(crate) struct SymmetricTensor6<T> {
    pub(crate) xx: T,
    pub(crate) yy: T,
    pub(crate) zz: T,
    pub(crate) xy: T,
    pub(crate) xz: T,
    pub(crate) yz: T,
}

#[inline]
pub(crate) fn linear_index(nx: usize, ny: usize, i: usize, j: usize, k: usize) -> usize {
    k * nx * ny + j * nx + i
}

#[inline]
pub(crate) fn derivative_x<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    i: usize,
    j: usize,
    k: usize,
    dx: T,
) -> Vector3<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let idx = linear_index(nx, ny, i, j, k);
    if nx == 1 {
        return Vector3::zeros();
    }
    if i == 0 {
        (velocity[idx + 1] - velocity[idx]) / dx
    } else if i + 1 == nx {
        (velocity[idx] - velocity[idx - 1]) / dx
    } else {
        let two = T::one() + T::one();
        (velocity[idx + 1] - velocity[idx - 1]) / (two * dx)
    }
}

#[inline]
pub(crate) fn derivative_y<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    i: usize,
    j: usize,
    k: usize,
    dy: T,
) -> Vector3<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let idx = linear_index(nx, ny, i, j, k);
    if ny == 1 {
        return Vector3::zeros();
    }
    if j == 0 {
        (velocity[idx + nx] - velocity[idx]) / dy
    } else if j + 1 == ny {
        (velocity[idx] - velocity[idx - nx]) / dy
    } else {
        let two = T::one() + T::one();
        (velocity[idx + nx] - velocity[idx - nx]) / (two * dy)
    }
}

#[inline]
pub(crate) fn derivative_z<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    nz: usize,
    i: usize,
    j: usize,
    k: usize,
    dz: T,
) -> Vector3<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let idx = linear_index(nx, ny, i, j, k);
    let plane = nx * ny;
    if nz == 1 {
        return Vector3::zeros();
    }
    if k == 0 {
        (velocity[idx + plane] - velocity[idx]) / dz
    } else if k + 1 == nz {
        (velocity[idx] - velocity[idx - plane]) / dz
    } else {
        let two = T::one() + T::one();
        (velocity[idx + plane] - velocity[idx - plane]) / (two * dz)
    }
}

#[inline]
pub(crate) fn velocity_gradient_tensor<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    nz: usize,
    i: usize,
    j: usize,
    k: usize,
    dx: T,
    dy: T,
    dz: T,
) -> [[T; 3]; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let du_dx = derivative_x(velocity, nx, ny, i, j, k, dx);
    let du_dy = derivative_y(velocity, nx, ny, i, j, k, dy);
    let du_dz = derivative_z(velocity, nx, ny, nz, i, j, k, dz);

    [
        [du_dx.x, du_dy.x, du_dz.x],
        [du_dx.y, du_dy.y, du_dz.y],
        [du_dx.z, du_dy.z, du_dz.z],
    ]
}

#[inline]
pub(crate) fn strain_components<T>(gradient: &[[T; 3]; 3]) -> SymmetricTensor6<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let half = T::one() / (T::one() + T::one());
    SymmetricTensor6 {
        xx: gradient[0][0],
        yy: gradient[1][1],
        zz: gradient[2][2],
        xy: (gradient[0][1] + gradient[1][0]) * half,
        xz: (gradient[0][2] + gradient[2][0]) * half,
        yz: (gradient[1][2] + gradient[2][1]) * half,
    }
}

#[inline]
pub(crate) fn symmetric_contract<T>(left: SymmetricTensor6<T>, right: SymmetricTensor6<T>) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    let two = T::one() + T::one();
    left.xx * right.xx
        + left.yy * right.yy
        + left.zz * right.zz
        + two * (left.xy * right.xy + left.xz * right.xz + left.yz * right.yz)
}

#[inline]
pub(crate) fn symmetric_matrix_contract<T>(left: SymmetricTensor6<T>, right: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    left.xx * right[0][0]
        + left.yy * right[1][1]
        + left.zz * right[2][2]
        + left.xy * (right[0][1] + right[1][0])
        + left.xz * (right[0][2] + right[2][0])
        + left.yz * (right[1][2] + right[2][1])
}

#[inline]
pub(crate) fn strain_magnitude<T>(gradient: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive + num_traits::Float,
{
    let strain = strain_components(gradient);
    let two = T::one() + T::one();
    num_traits::Float::sqrt(two * symmetric_contract(strain, strain))
}

#[inline]
pub(crate) fn matrix_contract<T>(left: &[[T; 3]; 3], right: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    let mut total = T::zero();
    for i in 0..3 {
        for j in 0..3 {
            total += left[i][j] * right[i][j];
        }
    }
    total
}

#[inline]
pub(crate) fn matrix_frobenius_norm_sq<T>(matrix: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    matrix_contract(matrix, matrix)
}

#[inline]
pub(crate) fn gradient_magnitude<T>(gradient: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive + num_traits::Float,
{
    let two = T::one() + T::one();
    num_traits::Float::sqrt(two * matrix_frobenius_norm_sq(gradient))
}

#[inline]
pub(crate) fn vorticity_magnitude<T>(gradient: &[[T; 3]; 3]) -> T
where
    T: RealField + Copy + FromPrimitive + num_traits::Float,
{
    let omega_x = gradient[2][1] - gradient[1][2];
    let omega_y = gradient[0][2] - gradient[2][0];
    let omega_z = gradient[1][0] - gradient[0][1];
    num_traits::Float::sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z)
}

#[inline]
pub(crate) fn matrix_square<T>(matrix: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let mut square = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut value = T::zero();
            for k in 0..3 {
                value += matrix[i][k] * matrix[k][j];
            }
            square[i][j] = value;
        }
    }
    square
}

#[inline]
pub(crate) fn symmetric_trace_free_part<T>(matrix: &[[T; 3]; 3]) -> SymmetricTensor6<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let half = T::one() / (T::one() + T::one());
    let third = T::one() / (T::one() + T::one() + T::one());
    let trace = matrix[0][0] + matrix[1][1] + matrix[2][2];

    SymmetricTensor6 {
        xx: matrix[0][0] - third * trace,
        yy: matrix[1][1] - third * trace,
        zz: matrix[2][2] - third * trace,
        xy: (matrix[0][1] + matrix[1][0]) * half,
        xz: (matrix[0][2] + matrix[2][0]) * half,
        yz: (matrix[1][2] + matrix[2][1]) * half,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_contract_counts_off_diagonals_twice() {
        let left = SymmetricTensor6 {
            xx: 0.0_f64,
            yy: 0.0,
            zz: 0.0,
            xy: 1.0,
            xz: 0.0,
            yz: 0.0,
        };
        let right = left;

        assert_eq!(symmetric_contract(left, right), 2.0);
    }

    #[test]
    fn symmetric_contract_matches_diagonal_norm() {
        let left = SymmetricTensor6 {
            xx: 3.0_f64,
            yy: 0.0,
            zz: 0.0,
            xy: 0.0,
            xz: 0.0,
            yz: 0.0,
        };
        let right = left;

        assert_eq!(symmetric_contract(left, right), 9.0);
    }

    #[test]
    fn symmetric_matrix_contract_matches_symmetric_contract() {
        let tensor = SymmetricTensor6 {
            xx: 1.0_f64,
            yy: 2.0,
            zz: 3.0,
            xy: 4.0,
            xz: 5.0,
            yz: 6.0,
        };
        let matrix = [[1.0_f64, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]];

        assert_eq!(
            symmetric_matrix_contract(tensor, &matrix),
            symmetric_contract(tensor, tensor)
        );
    }

    #[test]
    fn gradient_magnitude_matches_manual_frobenius_norm() {
        let gradient = [[1.0_f64, -2.0, 3.0], [4.0, 5.0, -6.0], [-7.0, 8.0, 9.0]];
        let expected = (2.0_f64 * 285.0).sqrt();

        assert_eq!(matrix_frobenius_norm_sq(&gradient), 285.0);
        assert_eq!(gradient_magnitude(&gradient), expected);
    }
}
