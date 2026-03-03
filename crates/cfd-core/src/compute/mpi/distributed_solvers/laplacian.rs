use super::traits::DistributedLinearOperator;
use super::vector::DistributedVector;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::decomposition::{DomainDecomposition, LocalSubdomain};
use crate::compute::mpi::error::MpiResult;
use crate::compute::mpi::ghost_cells::GhostCellManager;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// 2D Laplacian operator for distributed domains
pub struct DistributedLaplacian2D<T: RealField> {
    communicator: MpiCommunicator,
    subdomain: LocalSubdomain,
    ghost_manager: GhostCellManager<T>,
    dx: T,
    dy: T,
    nx_global: usize,
    ny_global: usize,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedLaplacian2D<T> {
    /// Create new distributed 2D Laplacian
    pub fn new(
        decomp: &DomainDecomposition,
        communicator: &MpiCommunicator,
        dx: T,
        dy: T,
    ) -> MpiResult<Self> {
        let ghost_manager = GhostCellManager::new(
            communicator.clone(),
            decomp.neighbors().clone(),
            1, // ghost layers
        )?;

        Ok(Self {
            communicator: communicator.clone(),
            subdomain: decomp.local_subdomain().clone(),
            ghost_manager,
            dx,
            dy,
            nx_global: decomp.global_extents().nx_global,
            ny_global: decomp.global_extents().ny_global,
        })
    }
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedLinearOperator<T>
    for DistributedLaplacian2D<T>
{
    fn apply(&self, x: &DistributedVector<T>, y: &mut DistributedVector<T>) -> MpiResult<()> {
        // Apply local Laplacian operator
        self.apply_local_laplacian(&x.local_data, &mut y.local_data);

        // Update ghost cells in the result
        if let Some(ref ghost_mgr) = y.ghost_manager {
            // Convert data layout for ghost cell exchange
            let mut velocity_u = vec![
                vec![nalgebra::Vector2::zeros(); self.subdomain.total_ny()];
                self.subdomain.total_nx()
            ];
            let mut velocity_v = vec![
                vec![nalgebra::Vector2::zeros(); self.subdomain.total_ny()];
                self.subdomain.total_nx()
            ];
            let mut pressure =
                vec![vec![T::zero(); self.subdomain.total_ny()]; self.subdomain.total_nx()];

            // Pack scalar field into pressure array for ghost exchange
            for i in 0..self.subdomain.nx_local {
                for j in 0..self.subdomain.ny_local {
                    let idx = i * self.subdomain.ny_local + j;
                    pressure[i + self.subdomain.ghost_layers][j + self.subdomain.ghost_layers] =
                        y.local_data[idx];
                }
            }

            ghost_mgr.update_ghost_cells(
                &mut velocity_u,
                &mut velocity_v,
                &mut pressure,
                &self.subdomain,
            )?;

            // Unpack ghost cells back to distributed vector
            for i in 0..self.subdomain.total_nx() {
                for j in 0..self.subdomain.total_ny() {
                    if !self.subdomain.is_owned_cell(i, j, 0) {
                        // This is a ghost cell - update the corresponding local index
                        let global_i = self.subdomain.local_to_global(i, j, 0).0;
                        let global_j = self.subdomain.local_to_global(i, j, 0).1;

                        if self.subdomain.owns_global_cell(global_i, global_j, 0) {
                            let local_i = global_i - self.subdomain.i_start_global;
                            let local_j = global_j - self.subdomain.j_start_global;
                            if local_i < self.subdomain.nx_local
                                && local_j < self.subdomain.ny_local
                            {
                                let idx = local_i * self.subdomain.ny_local + local_j;
                                y.local_data[idx] = pressure[i][j];
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn local_dimension(&self) -> usize {
        self.subdomain.nx_local * self.subdomain.ny_local
    }

    fn global_dimension(&self) -> usize {
        self.nx_global * self.ny_global
    }

    fn extract_diagonal(&self) -> DVector<T> {
        let nx = self.subdomain.nx_local;
        let ny = self.subdomain.ny_local;
        let dx_sq = self.dx * self.dx;
        let dy_sq = self.dy * self.dy;
        let mut diag = DVector::zeros(nx * ny);

        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let mut val = T::zero();

                // x-direction
                if i > 0 { val -= T::one() / dx_sq; }
                if i < nx - 1 { val -= T::one() / dx_sq; }

                // y-direction
                if j > 0 { val -= T::one() / dy_sq; }
                if j < ny - 1 { val -= T::one() / dy_sq; }

                diag[idx] = val;
            }
        }
        diag
    }

    fn assemble_local_matrix(&self) -> Option<nalgebra::DMatrix<T>> {
        let nx = self.subdomain.nx_local;
        let ny = self.subdomain.ny_local;
        let n = nx * ny;
        let dx_sq = self.dx * self.dx;
        let dy_sq = self.dy * self.dy;
        let mut mat = nalgebra::DMatrix::zeros(n, n);

        for i in 0..nx {
            for j in 0..ny {
                let row = i * ny + j;

                // Diagonal (calculated same as extract_diagonal)
                let mut diag_val = T::zero();

                // x-direction
                if i > 0 {
                    diag_val -= T::one() / dx_sq;
                    mat[(row, (i - 1) * ny + j)] = T::one() / dx_sq;
                }
                if i < nx - 1 {
                    diag_val -= T::one() / dx_sq;
                    mat[(row, (i + 1) * ny + j)] = T::one() / dx_sq;
                }

                // y-direction
                if j > 0 {
                    diag_val -= T::one() / dy_sq;
                    mat[(row, i * ny + (j - 1))] = T::one() / dy_sq;
                }
                if j < ny - 1 {
                    diag_val -= T::one() / dy_sq;
                    mat[(row, i * ny + (j + 1))] = T::one() / dy_sq;
                }

                mat[(row, row)] = diag_val;
            }
        }
        Some(mat)
    }
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedLaplacian2D<T> {
    /// Apply local Laplacian operator (interior points only)
    fn apply_local_laplacian(&self, x: &DVector<T>, y: &mut DVector<T>) {
        let nx = self.subdomain.nx_local;
        let ny = self.subdomain.ny_local;
        let dx_sq = self.dx * self.dx;
        let dy_sq = self.dy * self.dy;

        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let mut laplacian_val = T::zero();

                // x-direction second derivatives
                if i > 0 {
                    laplacian_val += (x[(i - 1) * ny + j] - x[idx]) / dx_sq;
                }
                if i < nx - 1 {
                    laplacian_val += (x[(i + 1) * ny + j] - x[idx]) / dx_sq;
                }

                // y-direction second derivatives
                if j > 0 {
                    laplacian_val += (x[i * ny + (j - 1)] - x[idx]) / dy_sq;
                }
                if j < ny - 1 {
                    laplacian_val += (x[i * ny + (j + 1)] - x[idx]) / dy_sq;
                }

                y[idx] = laplacian_val;
            }
        }
    }
}
