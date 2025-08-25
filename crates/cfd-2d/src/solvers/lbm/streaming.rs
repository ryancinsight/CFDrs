//! Streaming operators for LBM.
//!
//! This module handles the propagation of distribution functions
//! along lattice links.

use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;

/// Streaming operator for propagating distributions
pub struct StreamingOperator;

impl StreamingOperator {
    /// Perform streaming step (pull scheme)
    pub fn stream<T: RealField + Copy>(f_src: &Vec<Vec<[T; 9]>>, f_dst: &mut Vec<Vec<[T; 9]>>) {
        let ny = f_src.len();
        let nx = if ny > 0 { f_src[0].len() } else { 0 };

        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];

                    // Source indices with periodic boundary conditions
                    let src_i = ((i as i32 - ex + nx as i32) % nx as i32) as usize;
                    let src_j = ((j as i32 - ey + ny as i32) % ny as i32) as usize;

                    f_dst[j][i][q] = f_src[src_j][src_i][q];
                }
            }
        }
    }

    /// Perform streaming with non-periodic boundaries
    pub fn stream_with_boundaries<T: RealField + Copy>(
        f_src: &Vec<Vec<[T; 9]>>,
        f_dst: &mut Vec<Vec<[T; 9]>>,
        boundary_mask: &Vec<Vec<bool>>,
    ) {
        let ny = f_src.len();
        let nx = if ny > 0 { f_src[0].len() } else { 0 };

        for j in 0..ny {
            for i in 0..nx {
                if boundary_mask[j][i] {
                    // Skip boundary nodes (handled separately)
                    continue;
                }

                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];

                    let src_i = i as i32 - ex;
                    let src_j = j as i32 - ey;

                    // Check if source is within domain
                    if src_i >= 0 && src_i < nx as i32 && src_j >= 0 && src_j < ny as i32 {
                        let src_i = src_i as usize;
                        let src_j = src_j as usize;

                        if !boundary_mask[src_j][src_i] {
                            f_dst[j][i][q] = f_src[src_j][src_i][q];
                        }
                    }
                }
            }
        }
    }

    /// Push scheme streaming (alternative implementation)
    pub fn stream_push<T: RealField + Copy>(
        f_src: &Vec<Vec<[T; 9]>>,
        f_dst: &mut Vec<Vec<[T; 9]>>,
    ) {
        let ny = f_src.len();
        let nx = if ny > 0 { f_src[0].len() } else { 0 };

        // Initialize destination with zeros
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    f_dst[j][i][q] = T::zero();
                }
            }
        }

        // Push distributions to neighbors
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];

                    // Destination indices with periodic boundary conditions
                    let dst_i = ((i as i32 + ex + nx as i32) % nx as i32) as usize;
                    let dst_j = ((j as i32 + ey + ny as i32) % ny as i32) as usize;

                    f_dst[dst_j][dst_i][q] = f_src[j][i][q];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_preserves_mass() {
        let nx = 10;
        let ny = 10;

        // Initialize with uniform distribution
        let mut f_src = vec![vec![[1.0_f64 / 9.0; 9]; nx]; ny];
        let mut f_dst = vec![vec![[0.0_f64; 9]; nx]; ny];

        // Compute initial mass
        let initial_mass: f64 = f_src
            .iter()
            .flat_map(|row| row.iter())
            .flat_map(|cell| cell.iter())
            .sum();

        // Perform streaming
        StreamingOperator::stream(&f_src, &mut f_dst);

        // Compute final mass
        let final_mass: f64 = f_dst
            .iter()
            .flat_map(|row| row.iter())
            .flat_map(|cell| cell.iter())
            .sum();

        // Mass should be conserved
        assert!((initial_mass - final_mass).abs() < 1e-10);
    }
}
