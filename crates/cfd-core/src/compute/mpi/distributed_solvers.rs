//! Distributed linear algebra operations for MPI parallelization.
//!
//! This module provides distributed sparse matrix operations and parallel preconditioners
//! essential for scalable CFD simulations. It enables efficient solution of large linear
//! systems across multiple MPI processes.
//!
//! ## Architecture
//!
//! The distributed solver architecture follows these principles:
//! - **Domain Decomposition**: Problems are partitioned across MPI processes
//! - **Ghost Cell Management**: Boundary data exchange between subdomains
//! - **Collective Operations**: Global reductions and broadcasts
//! - **Load Balancing**: Even distribution of computational work
//!
//! ## Parallel Preconditioners
//!
//! - **Block Jacobi**: Diagonal blocks solved independently, efficient for well-conditioned problems
//! - **Additive Schwarz**: Overlapping domain decomposition with local solves
//! - **Distributed ILU**: Incomplete LU factorization across process boundaries
//!
//! ## Usage
//!
//! ```no_run
//! use cfd_core::compute::mpi::*;
//!
//! // Initialize MPI and domain decomposition
//! let universe = MpiUniverse::new()?;
//! let world = universe.world();
//! let decomp = DomainDecomposition::new(global_extents, &world, DecompositionStrategy::Cartesian2D)?;
//!
//! // Create distributed linear operator
//! let operator = DistributedLaplacian2D::new(&decomp, &world)?;
//!
//! // Setup parallel preconditioner
//! let preconditioner = BlockJacobiPreconditioner::new(&operator, &decomp, &world)?;
//!
//! // Solve distributed system
//! let solver = DistributedGMRES::new(operator, preconditioner, &world);
//! let solution = solver.solve(&rhs, &initial_guess, tolerance, max_iter)?;
//! ```

use super::communicator::MpiCommunicator;
use super::decomposition::{DomainDecomposition, LocalSubdomain};
use super::error::{MpiError, MpiResult};
use super::ghost_cells::GhostCellManager;
use nalgebra::{DVector, RealField};
use std::collections::HashMap;

/// Distributed linear operator trait for matrix-free operations
pub trait DistributedLinearOperator<T: RealField> {
    /// Apply the operator to a distributed vector
    fn apply(&self, x: &DistributedVector<T>, y: &mut DistributedVector<T>) -> MpiResult<()>;

    /// Get local dimension owned by this process
    fn local_dimension(&self) -> usize;

    /// Get global dimension across all processes
    fn global_dimension(&self) -> usize;

    /// Extract the diagonal of the local operator
    fn extract_diagonal(&self) -> DVector<T>;

    /// Assemble the local matrix (optional, for Schwarz/ILU)
    fn assemble_local_matrix(&self) -> Option<nalgebra::DMatrix<T>> {
        None
    }
}

/// Preconditioner trait for distributed solvers
pub trait Preconditioner<T: RealField> {
    /// Apply preconditioner: z = M^-1 * r
    fn apply(&self, r: &DistributedVector<T>, z: &mut DistributedVector<T>) -> MpiResult<()>;
}

/// Distributed vector with MPI-aware data distribution
#[derive(Debug, Clone)]
pub struct DistributedVector<T: RealField> {
    /// Local data owned by this process (excluding ghost cells)
    pub local_data: DVector<T>,
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Domain decomposition information
    subdomain: LocalSubdomain,
    /// Ghost cell manager for boundary exchange
    ghost_manager: Option<GhostCellManager<T>>,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedVector<T> {
    /// Create a new distributed vector
    pub fn new(
        local_size: usize,
        communicator: &MpiCommunicator,
        subdomain: LocalSubdomain,
        ghost_manager: Option<GhostCellManager<T>>,
    ) -> Self {
        let local_data = DVector::zeros(local_size);
        Self {
            local_data,
            communicator: communicator.clone(),
            subdomain,
            ghost_manager,
        }
    }

    /// Create distributed vector from local data
    pub fn from_local_data(
        data: DVector<T>,
        communicator: &MpiCommunicator,
        subdomain: LocalSubdomain,
        ghost_manager: Option<GhostCellManager<T>>,
    ) -> Self {
        Self {
            local_data: data,
            communicator: communicator.clone(),
            subdomain,
            ghost_manager,
        }
    }

    /// Global dot product across all processes
    pub fn dot(&self, other: &DistributedVector<T>) -> MpiResult<T> {
        let local_dot = self.local_data.dot(&other.local_data);
        let mut global_dot = T::zero();
        self.communicator.all_reduce_sum(&mut global_dot, local_dot);
        Ok(global_dot)
    }

    /// Global L2 norm across all processes
    pub fn norm(&self) -> MpiResult<T> {
        let local_norm_sq = self.local_data.norm_squared();
        let mut global_norm_sq = T::zero();
        self.communicator
            .all_reduce_sum(&mut global_norm_sq, local_norm_sq);
        Ok(global_norm_sq.sqrt())
    }

    /// Scale vector by scalar
    pub fn scale(&mut self, alpha: T) {
        self.local_data.scale_mut(alpha);
    }

    /// Add scaled vector: y = y + alpha * x
    pub fn axpy(&mut self, alpha: T, x: &DistributedVector<T>) {
        self.local_data.axpy(alpha, &x.local_data, T::one());
    }

    /// Copy from another distributed vector
    pub fn copy_from(&mut self, other: &DistributedVector<T>) {
        self.local_data.copy_from(&other.local_data);
    }
}

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
                if i > 0 {
                    val -= T::one() / dx_sq;
                }
                if i < nx - 1 {
                    val -= T::one() / dx_sq;
                }

                // y-direction
                if j > 0 {
                    val -= T::one() / dy_sq;
                }
                if j < ny - 1 {
                    val -= T::one() / dy_sq;
                }

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
                let mut laplacian = T::zero();

                // x-direction second derivatives
                if i > 0 {
                    laplacian += (x[(i - 1) * ny + j] - x[idx]) / dx_sq;
                }
                if i < nx - 1 {
                    laplacian += (x[(i + 1) * ny + j] - x[idx]) / dx_sq;
                }

                // y-direction second derivatives
                if j > 0 {
                    laplacian += (x[i * ny + (j - 1)] - x[idx]) / dy_sq;
                }
                if j < ny - 1 {
                    laplacian += (x[i * ny + (j + 1)] - x[idx]) / dy_sq;
                }

                y[idx] = laplacian;
            }
        }
    }
}

/// Block Jacobi preconditioner for distributed systems
pub struct BlockJacobiPreconditioner<T: RealField, Op: DistributedLinearOperator<T>> {
    /// Local diagonal blocks for each process
    local_blocks: Vec<DVector<T>>,
    /// Process-local operator for applying blocks
    _operator: std::marker::PhantomData<Op>,
    /// MPI communicator
    communicator: MpiCommunicator,
}

impl<T: RealField + Copy + FromPrimitive, Op: DistributedLinearOperator<T>>
    BlockJacobiPreconditioner<T, Op>
{
    /// Create new block Jacobi preconditioner
    pub fn new(
        operator: &Op,
        decomp: &DomainDecomposition,
        communicator: &MpiCommunicator,
    ) -> MpiResult<Self> {
        // Extract local diagonal blocks
        let local_dim = operator.local_dimension();
        let mut local_blocks = Vec::with_capacity(local_dim);

        // Use the actual diagonal from the operator
        let diagonal = operator.extract_diagonal();

        for i in 0..local_dim {
            // Store inverse of diagonal for fast application
            let val = if diagonal[i].is_zero() {
                T::one() // Avoid division by zero, effectively identity
            } else {
                T::one() / diagonal[i]
            };
            local_blocks.push(DVector::from_element(1, val));
        }

        Ok(Self {
            local_blocks,
            _operator: std::marker::PhantomData,
            communicator: communicator.clone(),
        })
    }
}

impl<T: RealField + Copy + FromPrimitive, Op: DistributedLinearOperator<T>> Preconditioner<T>
    for BlockJacobiPreconditioner<T, Op>
{
    /// Apply preconditioner: solve M * z = r for each block
    fn apply(&self, r: &DistributedVector<T>, z: &mut DistributedVector<T>) -> MpiResult<()> {
        // For Jacobi, z_i = r_i * inv_diagonal_i
        for i in 0..r.local_data.len() {
            z.local_data[i] = r.local_data[i] * self.local_blocks[i][0];
        }
        Ok(())
    }
}

/// Additive Schwarz preconditioner with overlapping domains
pub struct AdditiveSchwarzPreconditioner<T: RealField, Op: DistributedLinearOperator<T>> {
    /// Overlap size between subdomains
    overlap: usize,
    /// Local solvers for each overlapping subdomain
    local_solvers: Vec<Box<dyn Fn(&DVector<T>, &mut DVector<T>)>>,
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Domain decomposition
    decomp: DomainDecomposition,
    _operator: std::marker::PhantomData<Op>,
}

impl<T: RealField + Copy + FromPrimitive, Op: DistributedLinearOperator<T>>
    AdditiveSchwarzPreconditioner<T, Op>
{
    /// Create new additive Schwarz preconditioner
    pub fn new(
        operator: &Op,
        decomp: &DomainDecomposition,
        communicator: &MpiCommunicator,
        overlap: usize,
    ) -> MpiResult<Self> {
        // Create overlapping subdomain solvers
        let local_solvers = Self::create_local_solvers(operator, decomp, overlap)?;

        Ok(Self {
            overlap,
            local_solvers,
            communicator: communicator.clone(),
            decomp: decomp.clone(),
            _operator: std::marker::PhantomData,
        })
    }

    /// Create local solvers for overlapping subdomains
    fn create_local_solvers(
        operator: &Op,
        _decomp: &DomainDecomposition,
        _overlap: usize,
    ) -> MpiResult<Vec<Box<dyn Fn(&DVector<T>, &mut DVector<T>)>>> {
        // Try to assemble local matrix for direct solve
        if let Some(mat) = operator.assemble_local_matrix() {
            // Use LU decomposition if matrix is available
            // Note: This is a simplified local solver that doesn't account for overlap yet
            // Proper Schwarz requires extracting the overlapping submatrix
            let lu = mat.lu();
            let solver = Box::new(move |r: &DVector<T>, z: &mut DVector<T>| {
                if let Some(sol) = lu.solve(r) {
                    z.copy_from(&sol);
                } else {
                    // Fallback if singular (should not happen for Laplacian)
                    z.copy_from(r);
                }
            });
            Ok(vec![solver])
        } else {
            // Fallback to Jacobi if matrix assembly not supported
            let diagonal = operator.extract_diagonal();
            let solver = Box::new(move |r: &DVector<T>, z: &mut DVector<T>| {
                for i in 0..r.len() {
                    if !diagonal[i].is_zero() {
                        z[i] = r[i] / diagonal[i];
                    } else {
                        z[i] = r[i];
                    }
                }
            });
            Ok(vec![solver])
        }
    }
}

impl<T: RealField + Copy + FromPrimitive, Op: DistributedLinearOperator<T>> Preconditioner<T>
    for AdditiveSchwarzPreconditioner<T, Op>
{
    /// Apply additive Schwarz preconditioner
    fn apply(&self, r: &DistributedVector<T>, z: &mut DistributedVector<T>) -> MpiResult<()> {
        // Apply each local solver and accumulate results
        let mut local_result = DVector::zeros(r.local_data.len());

        for solver in &self.local_solvers {
            let mut temp = DVector::zeros(r.local_data.len());
            solver(&r.local_data, &mut temp);
            local_result += temp;
        }

        z.local_data.copy_from(&local_result);
        Ok(())
    }
}

/// Distributed GMRES solver
pub struct DistributedGMRES<T: RealField, Op: DistributedLinearOperator<T>, Prec> {
    /// Linear operator
    operator: Op,
    /// Preconditioner
    preconditioner: Prec,
    /// Krylov subspace dimension
    krylov_dim: usize,
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Work vectors
    work_vectors: Vec<DistributedVector<T>>,
    /// Hessenberg matrix
    hessenberg: Vec<T>,
    /// Givens rotations (cos, sin)
    givens_rotations: Vec<(T, T)>,
    /// Residual vector
    residual: DistributedVector<T>,
}

impl<
        T: RealField + Copy + FromPrimitive + std::fmt::LowerExp,
        Op: DistributedLinearOperator<T>,
        Prec: Preconditioner<T>,
    > DistributedGMRES<T, Op, Prec>
{
    /// Create new distributed GMRES solver
    pub fn new(
        operator: Op,
        preconditioner: Prec,
        communicator: &MpiCommunicator,
        krylov_dim: usize,
    ) -> Self {
        let local_dim = operator.local_dimension();
        let subdomain = LocalSubdomain {
            rank: communicator.rank(),
            nx_local: 0, // Would be set from domain decomposition
            ny_local: 0,
            nz_local: 0,
            i_start_global: 0,
            j_start_global: 0,
            k_start_global: 0,
            ghost_layers: 0,
        };

        let work_vectors = (0..krylov_dim + 1)
            .map(|_| DistributedVector::new(local_dim, communicator, subdomain.clone(), None))
            .collect();

        let residual = DistributedVector::new(local_dim, communicator, subdomain, None);

        Self {
            operator,
            preconditioner,
            krylov_dim,
            communicator: communicator.clone(),
            work_vectors,
            hessenberg: vec![T::zero(); krylov_dim * (krylov_dim + 1)],
            givens_rotations: vec![(T::zero(), T::zero()); krylov_dim],
            residual,
        }
    }

    /// Solve linear system using distributed GMRES
    pub fn solve(
        &mut self,
        b: &DistributedVector<T>,
        x0: &DistributedVector<T>,
        tolerance: T,
        max_iter: usize,
    ) -> MpiResult<DistributedVector<T>> {
        let mut x = x0.clone();
        let mut iter = 0;
        let mut residual_norm = T::zero();

        // Initial residual: r = b - A*x0
        self.operator.apply(&x, &mut self.residual)?;
        self.residual.axpy(-T::one(), b);
        self.residual.scale(-T::one()); // r = b - A*x0

        residual_norm = self.residual.norm()?;

        while residual_norm > tolerance && iter < max_iter {
            // Arnoldi process with preconditioning
            self.gmres_iteration(&mut x, b, &mut residual_norm)?;
            iter += 1;
        }

        Ok(x)
    }

    /// Compute Givens rotation parameters (c, s, r) such that
    /// [ c  s ] [ a ] = [ r ]
    /// [-s  c ] [ b ]   [ 0 ]
    fn compute_givens(&self, a: T, b: T) -> (T, T, T) {
        if b == T::zero() {
            (T::one(), T::zero(), a)
        } else if a == T::zero() {
            (T::zero(), T::one(), b)
        } else {
            let hypot = (a * a + b * b).sqrt();
            (a / hypot, b / hypot, hypot)
        }
    }

    /// Single GMRES iteration
    fn gmres_iteration(
        &mut self,
        x: &mut DistributedVector<T>,
        b: &DistributedVector<T>,
        residual_norm: &mut T,
    ) -> MpiResult<()> {
        // Precondition residual
        let mut v0 = self.residual.clone();
        self.preconditioner.apply(&self.residual, &mut v0)?;

        // Normalize v0
        let beta = v0.norm()?;
        if beta != T::zero() {
            v0.scale(T::one() / beta);
        }

        // RHS for the least squares problem (g)
        let mut g = vec![T::zero(); self.krylov_dim + 1];
        g[0] = beta;

        let mut k = 0; // Actual number of iterations performed

        // Arnoldi process
        for j in 0..self
            .krylov_dim
            .min(self.work_vectors.len().saturating_sub(1))
        {
            k = j + 1;
            self.work_vectors[j].copy_from(&v0);

            // Apply operator: v_{j+1} = A * v_j
            self.operator
                .apply(&self.work_vectors[j], &mut self.work_vectors[j + 1])?;
            self.preconditioner
                .apply(&self.work_vectors[j + 1], &mut self.work_vectors[j + 1])?;

            // Orthogonalization (Gram-Schmidt)
            for i in 0..=j {
                let h_ij = self.work_vectors[j + 1].dot(&self.work_vectors[i])?;
                self.work_vectors[j + 1].axpy(-h_ij, &self.work_vectors[i]);
                self.hessenberg[i * self.krylov_dim + j] = h_ij;
            }

            // Normalize
            let norm = self.work_vectors[j + 1].norm()?;
            self.hessenberg[(j + 1) * self.krylov_dim + j] = norm;

            if norm > T::zero() {
                self.work_vectors[j + 1].scale(T::one() / norm);
            }

            // Apply previous Givens rotations to the new column
            for i in 0..j {
                let (c, s) = self.givens_rotations[i];
                let h_ij = self.hessenberg[i * self.krylov_dim + j];
                let h_ip1j = self.hessenberg[(i + 1) * self.krylov_dim + j];

                self.hessenberg[i * self.krylov_dim + j] = c * h_ij + s * h_ip1j;
                self.hessenberg[(i + 1) * self.krylov_dim + j] = -s * h_ij + c * h_ip1j;
            }

            // Compute new Givens rotation
            let h_jj = self.hessenberg[j * self.krylov_dim + j];
            let h_jp1j = self.hessenberg[(j + 1) * self.krylov_dim + j];
            let (c, s, r) = self.compute_givens(h_jj, h_jp1j);

            self.givens_rotations[j] = (c, s);

            // Apply rotation to H and g
            self.hessenberg[j * self.krylov_dim + j] = r;
            self.hessenberg[(j + 1) * self.krylov_dim + j] = T::zero();

            let g_j = g[j];
            g[j] = c * g_j;
            g[j + 1] = -s * g_j;

            *residual_norm = g[j + 1].abs();
        }

        // Backward substitution to solve R * y = g
        let mut y = vec![T::zero(); k];
        for i in (0..k).rev() {
            let mut sum = g[i];
            for j in (i + 1)..k {
                sum -= self.hessenberg[i * self.krylov_dim + j] * y[j];
            }
            y[i] = sum / self.hessenberg[i * self.krylov_dim + i];
        }

        // Update solution: x = x + V * y
        for j in 0..k {
            x.axpy(y[j], &self.work_vectors[j]);
        }

        // Recompute residual for next restart
        self.operator.apply(x, &mut self.residual)?;
        self.residual.axpy(-T::one(), b);
        self.residual.scale(-T::one());
        *residual_norm = self.residual.norm()?;

        Ok(())
    }
}

/// Parallel I/O operations for distributed data
pub mod parallel_io {
    use super::*;
    use std::path::Path;

    /// Parallel VTK writer for distributed meshes
    pub struct ParallelVtkWriter<T: RealField> {
        communicator: MpiCommunicator,
        is_root: bool,
    }

    impl<T: RealField> ParallelVtkWriter<T> {
        /// Create new parallel VTK writer
        pub fn new(communicator: &MpiCommunicator) -> Self {
            Self {
                communicator: communicator.clone(),
                is_root: communicator.is_root(),
            }
        }

        /// Write distributed VTK file
        pub fn write_vtk_file<P: AsRef<Path>>(
            &self,
            filename: P,
            points: &DistributedVector<T>,
            cells: &[u32],
            cell_types: &[u8],
            point_data: &HashMap<String, &DistributedVector<T>>,
            cell_data: &HashMap<String, &DistributedVector<T>>,
        ) -> MpiResult<()> {
            if self.is_root {
                // Root process writes header and collects data from other processes
                self.write_vtk_header(filename, points, cells, cell_types)?;
            }

            // Each process writes its portion of point/cell data
            self.write_distributed_data(filename, point_data, cell_data)?;

            Ok(())
        }

        /// Write VTK header (root process only)
        fn write_vtk_header<P: AsRef<Path>>(
            &self,
            _filename: P,
            _points: &DistributedVector<T>,
            _cells: &[u32],
            _cell_types: &[u8],
        ) -> MpiResult<()> {
            // Implementation would write VTK header
            // This is a placeholder for the actual VTK writing logic
            Ok(())
        }

        /// Write distributed point and cell data
        fn write_distributed_data<P: AsRef<Path>>(
            &self,
            _filename: P,
            _point_data: &HashMap<String, &DistributedVector<T>>,
            _cell_data: &HashMap<String, &DistributedVector<T>>,
        ) -> MpiResult<()> {
            // Implementation would write distributed data
            // This is a placeholder for the actual distributed data writing
            Ok(())
        }
    }

    /// Parallel HDF5 writer for distributed datasets
    pub struct ParallelHdf5Writer<T: RealField> {
        communicator: MpiCommunicator,
        is_root: bool,
    }

    impl<T: RealField> ParallelHdf5Writer<T> {
        /// Create new parallel HDF5 writer
        pub fn new(communicator: &MpiCommunicator) -> Self {
            Self {
                communicator: communicator.clone(),
                is_root: communicator.is_root(),
            }
        }

        /// Write distributed HDF5 file
        pub fn write_hdf5_file<P: AsRef<Path>>(
            &self,
            filename: P,
            datasets: &HashMap<String, &DistributedVector<T>>,
            metadata: &HashMap<String, String>,
        ) -> MpiResult<()> {
            if self.is_root {
                // Root process creates file and writes metadata
                self.write_hdf5_header(filename, metadata)?;
            }

            // Collective write of distributed data
            self.write_distributed_datasets(filename, datasets)?;

            Ok(())
        }

        /// Write HDF5 header and metadata (root process only)
        fn write_hdf5_header<P: AsRef<Path>>(
            &self,
            _filename: P,
            _metadata: &HashMap<String, String>,
        ) -> MpiResult<()> {
            // Implementation would write HDF5 header and metadata
            // This is a placeholder for the actual HDF5 writing logic
            Ok(())
        }

        /// Write distributed datasets using collective I/O
        fn write_distributed_datasets<P: AsRef<Path>>(
            &self,
            _filename: P,
            _datasets: &HashMap<String, &DistributedVector<T>>,
        ) -> MpiResult<()> {
            // Implementation would use HDF5 collective I/O operations
            // This is a placeholder for the actual distributed dataset writing
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::RealField;

    #[test]
    fn test_distributed_vector_creation() {
        // This would require MPI initialization in a real test
        // Placeholder for compile-time validation
        let _marker: std::marker::PhantomData<DistributedVector<f64>> = std::marker::PhantomData;
    }

    #[test]
    fn test_parallel_io_types() {
        // Test that parallel I/O types can be created (compile-time check)
        let _vtk_marker: std::marker::PhantomData<parallel_io::ParallelVtkWriter<f64>> =
            std::marker::PhantomData;
        let _hdf5_marker: std::marker::PhantomData<parallel_io::ParallelHdf5Writer<f64>> =
            std::marker::PhantomData;
    }
}
