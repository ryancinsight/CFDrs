//! f64 SWAR operations

use super::SwarOps;
use crate::error::Result;

impl SwarOps {
    /// Add f64 arrays with loop unrolling
    #[inline]
    pub fn add_f64_arrays(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if self.unroll_factor == 4 {
            let len = a.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (res_main, res_rem) = result.split_at_mut(split_idx);
            let (a_main, a_rem) = a.split_at(split_idx);
            let (b_main, b_rem) = b.split_at(split_idx);

            for ((ca, cb), cr) in a_main
                .chunks_exact(4)
                .zip(b_main.chunks_exact(4))
                .zip(res_main.chunks_exact_mut(4))
            {
                cr[0] = ca[0] + cb[0];
                cr[1] = ca[1] + cb[1];
                cr[2] = ca[2] + cb[2];
                cr[3] = ca[3] + cb[3];
            }

            for ((x, y), z) in a_rem.iter().zip(b_rem.iter()).zip(res_rem.iter_mut()) {
                *z = *x + *y;
            }
        } else {
            // Fallback for non-standard unroll factors
            for ((x, y), z) in a.iter().zip(b.iter()).zip(result.iter_mut()) {
                *z = *x + *y;
            }
        }

        Ok(())
    }

    /// Subtract f64 arrays with loop unrolling
    #[inline]
    pub fn sub_f64_arrays(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if self.unroll_factor == 4 {
            let len = a.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (res_main, res_rem) = result.split_at_mut(split_idx);
            let (a_main, a_rem) = a.split_at(split_idx);
            let (b_main, b_rem) = b.split_at(split_idx);

            for ((ca, cb), cr) in a_main
                .chunks_exact(4)
                .zip(b_main.chunks_exact(4))
                .zip(res_main.chunks_exact_mut(4))
            {
                cr[0] = ca[0] - cb[0];
                cr[1] = ca[1] - cb[1];
                cr[2] = ca[2] - cb[2];
                cr[3] = ca[3] - cb[3];
            }

            for ((x, y), z) in a_rem.iter().zip(b_rem.iter()).zip(res_rem.iter_mut()) {
                *z = *x - *y;
            }
        } else {
            for ((x, y), z) in a.iter().zip(b.iter()).zip(result.iter_mut()) {
                *z = *x - *y;
            }
        }

        Ok(())
    }

    /// Multiply f64 arrays with loop unrolling
    #[inline]
    pub fn mul_f64_arrays(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if self.unroll_factor == 4 {
            let len = a.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (res_main, res_rem) = result.split_at_mut(split_idx);
            let (a_main, a_rem) = a.split_at(split_idx);
            let (b_main, b_rem) = b.split_at(split_idx);

            for ((ca, cb), cr) in a_main
                .chunks_exact(4)
                .zip(b_main.chunks_exact(4))
                .zip(res_main.chunks_exact_mut(4))
            {
                cr[0] = ca[0] * cb[0];
                cr[1] = ca[1] * cb[1];
                cr[2] = ca[2] * cb[2];
                cr[3] = ca[3] * cb[3];
            }

            for ((x, y), z) in a_rem.iter().zip(b_rem.iter()).zip(res_rem.iter_mut()) {
                *z = *x * *y;
            }
        } else {
            for ((x, y), z) in a.iter().zip(b.iter()).zip(result.iter_mut()) {
                *z = *x * *y;
            }
        }

        Ok(())
    }

    /// Divide f64 arrays
    #[inline]
    pub fn div_f64_arrays(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        for i in 0..a.len() {
            result[i] = a[i] / b[i];
        }
        Ok(())
    }

    /// Fused multiply-add for f64
    #[inline]
    pub fn fma_f64_arrays(
        &self,
        a: &[f64],
        b: &[f64],
        c: &[f64],
        result: &mut [f64],
    ) -> Result<()> {
        if self.unroll_factor == 4 {
            let len = a.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (res_main, res_rem) = result.split_at_mut(split_idx);
            let (a_main, a_rem) = a.split_at(split_idx);
            let (b_main, b_rem) = b.split_at(split_idx);
            let (c_main, c_rem) = c.split_at(split_idx);

            for (((ca, cb), cc), cr) in a_main
                .chunks_exact(4)
                .zip(b_main.chunks_exact(4))
                .zip(c_main.chunks_exact(4))
                .zip(res_main.chunks_exact_mut(4))
            {
                cr[0] = ca[0].mul_add(cb[0], cc[0]);
                cr[1] = ca[1].mul_add(cb[1], cc[1]);
                cr[2] = ca[2].mul_add(cb[2], cc[2]);
                cr[3] = ca[3].mul_add(cb[3], cc[3]);
            }

            for (((x, y), z), w) in a_rem
                .iter()
                .zip(b_rem.iter())
                .zip(c_rem.iter())
                .zip(res_rem.iter_mut())
            {
                *w = x.mul_add(*y, *z);
            }
        } else {
            for (((x, y), z), w) in a.iter().zip(b.iter()).zip(c.iter()).zip(result.iter_mut()) {
                *w = x.mul_add(*y, *z);
            }
        }

        Ok(())
    }

    /// Scale f64 array by scalar
    #[inline]
    pub fn scale_f64_array(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        if self.unroll_factor == 4 {
            let len = input.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (res_main, res_rem) = result.split_at_mut(split_idx);
            let (in_main, in_rem) = input.split_at(split_idx);

            for (inp, res) in in_main.chunks_exact(4).zip(res_main.chunks_exact_mut(4)) {
                res[0] = inp[0] * scalar;
                res[1] = inp[1] * scalar;
                res[2] = inp[2] * scalar;
                res[3] = inp[3] * scalar;
            }

            for (x, y) in in_rem.iter().zip(res_rem.iter_mut()) {
                *y = *x * scalar;
            }
        } else {
            for (x, y) in input.iter().zip(result.iter_mut()) {
                *y = *x * scalar;
            }
        }

        Ok(())
    }

    /// Dot product for f64
    #[inline]
    pub fn dot_f64_arrays(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        let mut sum = 0.0;
        if self.unroll_factor == 4 {
            let len = a.len();
            let chunk_count = len / 4;
            let split_idx = chunk_count * 4;

            let (a_main, a_rem) = a.split_at(split_idx);
            let (b_main, b_rem) = b.split_at(split_idx);

            for (ca, cb) in a_main.chunks_exact(4).zip(b_main.chunks_exact(4)) {
                sum += ca[0] * cb[0] + ca[1] * cb[1] + ca[2] * cb[2] + ca[3] * cb[3];
            }

            for (x, y) in a_rem.iter().zip(b_rem.iter()) {
                sum += *x * *y;
            }
        } else {
            for (x, y) in a.iter().zip(b.iter()) {
                sum += *x * *y;
            }
        }

        Ok(sum)
    }
}
