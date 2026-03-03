use crate::compute::mpi::error::{MpiError, MpiResult};
use crate::compute::mpi::{LocalSubdomain, Rank};

#[derive(Clone, Copy)]
pub(super) struct PartitionRegion {
    pub(super) i_start: usize,
    pub(super) j_start: usize,
    pub(super) k_start: usize,
    pub(super) nx: usize,
    pub(super) ny: usize,
    pub(super) nz: usize,
}

#[derive(Clone, Copy)]
pub(super) enum PartitionDimension {
    X,
    Y,
    Z,
}

impl PartitionRegion {
    pub(super) fn select_split_dimension(&self) -> Option<PartitionDimension> {
        let mut best = None;
        let mut best_len = 0usize;
        if self.nx > 1 && self.nx >= best_len {
            best = Some(PartitionDimension::X);
            best_len = self.nx;
        }
        if self.ny > 1 && self.ny > best_len {
            best = Some(PartitionDimension::Y);
            best_len = self.ny;
        }
        if self.nz > 1 && self.nz > best_len {
            best = Some(PartitionDimension::Z);
        }
        best
    }

    pub(super) fn split(
        &self,
        dim: PartitionDimension,
        left_weight: f64,
        right_weight: f64,
    ) -> MpiResult<(PartitionRegion, PartitionRegion)> {
        let total_weight = left_weight + right_weight;
        if total_weight <= 0.0 {
            return Err(MpiError::DecompositionError(
                "Partition weights must be positive".to_string(),
            ));
        }
        let (len, min_len) = match dim {
            PartitionDimension::X => (self.nx, 1usize),
            PartitionDimension::Y => (self.ny, 1usize),
            PartitionDimension::Z => (self.nz, 1usize),
        };
        if len <= 1 {
            return Err(MpiError::DecompositionError(
                "Cannot split dimension with length <= 1".to_string(),
            ));
        }
        let mut left_len =
            ((left_weight / total_weight) * len as f64).round() as usize;
        if left_len < min_len {
            left_len = min_len;
        }
        if left_len >= len {
            left_len = len - 1;
        }
        let right_len = len - left_len;
        let left = match dim {
            PartitionDimension::X => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: left_len,
                ny: self.ny,
                nz: self.nz,
            },
            PartitionDimension::Y => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: self.nx,
                ny: left_len,
                nz: self.nz,
            },
            PartitionDimension::Z => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: self.nx,
                ny: self.ny,
                nz: left_len,
            },
        };
        let right = match dim {
            PartitionDimension::X => PartitionRegion {
                i_start: self.i_start + left_len,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: right_len,
                ny: self.ny,
                nz: self.nz,
            },
            PartitionDimension::Y => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start + left_len,
                k_start: self.k_start,
                nx: self.nx,
                ny: right_len,
                nz: self.nz,
            },
            PartitionDimension::Z => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start + left_len,
                nx: self.nx,
                ny: self.ny,
                nz: right_len,
            },
        };
        Ok((left, right))
    }
}

/// Partition 1D domain by weights across ranks
pub(super) fn partition_1d_by_weights(
    nx: usize,
    weights: &[usize],
    size: crate::compute::mpi::Size,
) -> MpiResult<Vec<(usize, usize)>> {
    if weights.len() != nx {
        return Err(MpiError::DecompositionError(
            "Weight vector length must match nx".to_string(),
        ));
    }
    let total_weight: usize = weights.iter().sum();
    if total_weight == 0 {
        return Err(MpiError::DecompositionError(
            "Total weight must be positive".to_string(),
        ));
    }
    let mut partitions = Vec::with_capacity(size as usize);
    let mut start = 0usize;
    let mut acc_weight = 0usize;
    let mut target_weight = total_weight as f64 / size as f64;
    let mut current_rank = 0usize;
    for i in 0..nx {
        acc_weight += weights[i];
        let should_split = acc_weight as f64 >= target_weight && current_rank + 1 < size as usize;
        if should_split {
            let length = i + 1 - start;
            partitions.push((start, length));
            start = i + 1;
            current_rank += 1;
            target_weight =
                (total_weight as f64 * (current_rank as f64 + 1.0)) / size as f64;
        }
    }
    if start < nx {
        partitions.push((start, nx - start));
    }
    if partitions.len() != size as usize {
        return Err(MpiError::DecompositionError(
            "Failed to partition domain across ranks".to_string(),
        ));
    }
    Ok(partitions)
}

/// Recursive coordinate bisection partitioning
pub(super) fn partition_rcb(
    global: &crate::compute::mpi::GlobalExtents,
    weights: &[f64],
) -> MpiResult<Vec<LocalSubdomain>> {
    let size = weights.len();
    if size == 0 {
        return Err(MpiError::DecompositionError(
            "Weights must not be empty".to_string(),
        ));
    }
    if global.nx_global == 0 || global.ny_global == 0 || global.nz_global == 0 {
        return Err(MpiError::DecompositionError(
            "Global extents must be non-zero in all dimensions".to_string(),
        ));
    }
    let total_weight: f64 = weights.iter().sum();
    if total_weight <= 0.0 {
        return Err(MpiError::DecompositionError(
            "Total weight must be positive".to_string(),
        ));
    }
    let mut partitions = Vec::with_capacity(size);
    let ranks: Vec<Rank> = (0..size as i32).collect();
    let region = PartitionRegion {
        i_start: 0,
        j_start: 0,
        k_start: 0,
        nx: global.nx_global,
        ny: global.ny_global,
        nz: global.nz_global,
    };
    partition_rcb_recursive(region, &ranks, weights, &mut partitions)?;
    if partitions.len() != size {
        return Err(MpiError::DecompositionError(
            "Recursive bisection did not produce all partitions".to_string(),
        ));
    }
    Ok(partitions)
}

/// Recursive helper for RCB partitioning
fn partition_rcb_recursive(
    region: PartitionRegion,
    ranks: &[Rank],
    weights: &[f64],
    partitions: &mut Vec<LocalSubdomain>,
) -> MpiResult<()> {
    if ranks.len() != weights.len() {
        return Err(MpiError::DecompositionError(
            "Rank list and weights must match".to_string(),
        ));
    }
    if ranks.is_empty() {
        return Err(MpiError::DecompositionError(
            "Ranks must not be empty".to_string(),
        ));
    }
    if ranks.len() == 1 {
        partitions.push(LocalSubdomain {
            rank: ranks[0],
            nx_local: region.nx,
            ny_local: region.ny,
            nz_local: region.nz,
            i_start_global: region.i_start,
            j_start_global: region.j_start,
            k_start_global: region.k_start,
            ghost_layers: 1,
        });
        return Ok(());
    }
    let total_weight: f64 = weights.iter().sum();
    if total_weight <= 0.0 {
        return Err(MpiError::DecompositionError(
            "Total weight must be positive".to_string(),
        ));
    }
    let split_dim = region.select_split_dimension().ok_or_else(|| {
        MpiError::DecompositionError(
            "Cannot split region with no available dimension".to_string(),
        )
    })?;
    let mut cumulative = 0.0;
    let mut split_index = 0usize;
    for (i, weight) in weights.iter().enumerate() {
        cumulative += *weight;
        if cumulative >= total_weight / 2.0 {
            split_index = i + 1;
            break;
        }
    }
    if split_index == 0 || split_index >= ranks.len() {
        split_index = ranks.len() / 2;
    }
    let left_weights = &weights[..split_index];
    let right_weights = &weights[split_index..];
    let left_weight: f64 = left_weights.iter().sum();
    let right_weight: f64 = right_weights.iter().sum();
    if left_weight <= 0.0 || right_weight <= 0.0 {
        return Err(MpiError::DecompositionError(
            "Partition weights must be positive".to_string(),
        ));
    }

    let (left_region, right_region) =
        region.split(split_dim, left_weight, right_weight)?;

    partition_rcb_recursive(left_region, &ranks[..split_index], left_weights, partitions)?;
    partition_rcb_recursive(
        right_region,
        &ranks[split_index..],
        right_weights,
        partitions,
    )?;
    Ok(())
}
