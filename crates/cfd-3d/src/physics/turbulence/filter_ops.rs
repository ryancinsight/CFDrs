use nalgebra::Vector3;
use num_traits::FromPrimitive;

use super::field_ops::SymmetricTensor6;

#[derive(Clone, Copy, Debug)]
pub(crate) struct FilterMoments<T: cfd_mesh::domain::core::Scalar + nalgebra::RealField + Copy> {
    pub(crate) uu: T,
    pub(crate) vv: T,
    pub(crate) ww: T,
    pub(crate) uv: T,
    pub(crate) uw: T,
    pub(crate) vw: T,
}

impl<T: cfd_mesh::domain::core::Scalar + nalgebra::RealField + Copy> FilterMoments<T> {
    #[inline]
    pub(crate) fn zero() -> Self {
        Self {
            uu: T::zero(),
            vv: T::zero(),
            ww: T::zero(),
            uv: T::zero(),
            uw: T::zero(),
            vw: T::zero(),
        }
    }
}

#[inline]
fn visit_box_stencil<F>(nx: usize, ny: usize, nz: usize, i: usize, j: usize, k: usize, mut visit: F)
where
    F: FnMut(usize),
{
    let plane = nx * ny;

    for dk in [0isize, -1, 1] {
        let nk = k as isize + dk;
        if nk < 0 || nk >= nz as isize {
            continue;
        }
        for dj in [0isize, -1, 1] {
            let nj = j as isize + dj;
            if nj < 0 || nj >= ny as isize {
                continue;
            }
            let base = nk as usize * plane + nj as usize * nx;
            for di in [0isize, -1, 1] {
                let ni = i as isize + di;
                if ni < 0 || ni >= nx as isize {
                    continue;
                }
                visit(base + ni as usize);
            }
        }
    }
}

#[inline]
pub(crate) fn box_filter_velocity_at<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    nz: usize,
    i: usize,
    j: usize,
    k: usize,
) -> Vector3<T>
where
    T: cfd_mesh::domain::core::Scalar + nalgebra::RealField + Copy + FromPrimitive,
{
    let mut sum = Vector3::zeros();
    let mut count = 0usize;

    visit_box_stencil(nx, ny, nz, i, j, k, |idx| {
        sum += velocity[idx];
        count += 1;
    });

    let count_t =
        <T as FromPrimitive>::from_usize(count).expect("box filter stencil size is representable");
    sum / count_t
}

#[inline]
pub(crate) fn box_filter_moments_at<T>(
    velocity: &[Vector3<T>],
    nx: usize,
    ny: usize,
    nz: usize,
    i: usize,
    j: usize,
    k: usize,
) -> FilterMoments<T>
where
    T: cfd_mesh::domain::core::Scalar + nalgebra::RealField + Copy + FromPrimitive,
{
    let mut moments = FilterMoments::zero();
    let mut count = 0usize;

    visit_box_stencil(nx, ny, nz, i, j, k, |idx| {
        let v = velocity[idx];
        moments.uu += v.x * v.x;
        moments.vv += v.y * v.y;
        moments.ww += v.z * v.z;
        moments.uv += v.x * v.y;
        moments.uw += v.x * v.z;
        moments.vw += v.y * v.z;
        count += 1;
    });

    let count_t =
        <T as FromPrimitive>::from_usize(count).expect("box filter stencil size is representable");
    let inv = T::one() / count_t;
    moments.uu *= inv;
    moments.vv *= inv;
    moments.ww *= inv;
    moments.uv *= inv;
    moments.uw *= inv;
    moments.vw *= inv;
    moments
}

#[inline]
pub(crate) fn resolved_stress_tensor<T>(
    moments: FilterMoments<T>,
    filtered_velocity: Vector3<T>,
) -> SymmetricTensor6<T>
where
    T: cfd_mesh::domain::core::Scalar + nalgebra::RealField + Copy + FromPrimitive,
{
    SymmetricTensor6 {
        xx: moments.uu - filtered_velocity.x * filtered_velocity.x,
        yy: moments.vv - filtered_velocity.y * filtered_velocity.y,
        zz: moments.ww - filtered_velocity.z * filtered_velocity.z,
        xy: moments.uv - filtered_velocity.x * filtered_velocity.y,
        xz: moments.uw - filtered_velocity.x * filtered_velocity.z,
        yz: moments.vw - filtered_velocity.y * filtered_velocity.z,
    }
}
