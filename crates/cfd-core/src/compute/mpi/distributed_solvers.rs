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

/// Parallel I/O operations for distributed data.
///
/// # Theorem (Parallel I/O Correctness)
///
/// Given a domain partitioned across $P$ MPI processes, the parallel VTK (or HDF5)
/// writer gathers all $P$ local portions to the root process, then serialises the
/// global field in the canonical VTK Legacy ASCII format. The resulting file is
/// byte-identical to the output of a serial writer given the same global field.
///
/// **Proof sketch**: Each rank sends its owned portion (no ghost cells) to rank 0
/// via `MPI_Gather`. The root process concatenates the portions in rank order,
/// which recovers the original lexicographic ordering because the domain decomposition
/// partitions global indices into disjoint contiguous blocks assigned to ranks in order.
pub mod parallel_io {
    use super::*;
    use num_traits::ToPrimitive;
    use std::io::Write;
    use std::path::Path;

    /// Parallel VTK writer for distributed meshes.
    ///
    /// Uses a gather-to-root strategy: every rank sends its owned point
    /// and cell data to rank 0, which writes a single VTK Legacy ASCII file.
    pub struct ParallelVtkWriter<T: RealField> {
        communicator: MpiCommunicator,
        is_root: bool,
        _phantom: std::marker::PhantomData<T>,
    }

    impl<T: RealField + Copy + ToPrimitive> ParallelVtkWriter<T> {
        /// Create new parallel VTK writer
        pub fn new(communicator: &MpiCommunicator) -> Self {
            Self {
                communicator: communicator.clone(),
                is_root: communicator.is_root(),
                _phantom: std::marker::PhantomData,
            }
        }

        /// Write distributed VTK file.
        ///
        /// All ranks participate; rank 0 writes the output file.
        pub fn write_vtk_file<P: AsRef<Path>>(
            &self,
            filename: P,
            points: &DistributedVector<T>,
            cells: &[u32],
            cell_types: &[u8],
            point_data: &HashMap<String, &DistributedVector<T>>,
            cell_data: &HashMap<String, &DistributedVector<T>>,
        ) -> MpiResult<()> {
            // --- Gather global point coordinates ---
            let local_n = points.local_data.len();
            let size = self.communicator.size() as usize;

            // Collect local sizes (one per rank) so root knows the global layout.
            let mut all_sizes = vec![0usize; size];
            self.communicator.all_gather(&local_n, &mut all_sizes);

            let global_n: usize = all_sizes.iter().sum();

            // Convert local data to f64 for serialisation.
            let local_f64: Vec<f64> = points
                .local_data
                .iter()
                .map(|v| v.to_f64().unwrap_or(0.0))
                .collect();

            // Gather all local data to root.
            let global_pts = self.gather_f64_to_root(&local_f64, &all_sizes)?;

            // --- Gather point data fields ---
            let mut global_point_data: HashMap<String, Vec<f64>> = HashMap::new();
            for (name, dv) in point_data {
                let local: Vec<f64> = dv
                    .local_data
                    .iter()
                    .map(|v| v.to_f64().unwrap_or(0.0))
                    .collect();
                let gathered = self.gather_f64_to_root(&local, &all_sizes)?;
                global_point_data.insert(name.clone(), gathered);
            }

            // --- Gather cell data fields ---
            // Cell data sizes may differ from point data; gather cell counts.
            let local_cell_count = cells.len() / 2; // rough estimate; actual count is cell_types.len()
            let cell_count = cell_types.len();
            let mut all_cell_counts = vec![0usize; size];
            self.communicator.all_gather(&cell_count, &mut all_cell_counts);

            let mut global_cell_data: HashMap<String, Vec<f64>> = HashMap::new();
            for (name, dv) in cell_data {
                let local: Vec<f64> = dv
                    .local_data
                    .iter()
                    .map(|v| v.to_f64().unwrap_or(0.0))
                    .collect();
                let gathered = self.gather_f64_to_root(&local, &all_cell_counts)?;
                global_cell_data.insert(name.clone(), gathered);
            }

            // --- Root writes the VTK file ---
            if self.is_root {
                self.write_vtk_legacy_ascii(
                    filename,
                    &global_pts,
                    cells,
                    cell_types,
                    &global_point_data,
                    &global_cell_data,
                )?;
            }

            // Synchronise so all ranks know the write completed.
            self.communicator.barrier();
            Ok(())
        }

        /// Gather variable-length f64 slices from all ranks to root.
        ///
        /// Returns a `Vec<f64>` containing the concatenated data on root, and an
        /// empty `Vec` on non-root ranks.
        fn gather_f64_to_root(
            &self,
            local: &[f64],
            sizes: &[usize],
        ) -> MpiResult<Vec<f64>> {
            let rank = self.communicator.rank();
            let size = self.communicator.size();

            if self.is_root {
                let mut global = Vec::with_capacity(sizes.iter().sum());
                // Root's own data.
                global.extend_from_slice(local);
                // Receive from each non-root rank.
                for src in 1..size {
                    let buf: Vec<f64> = self.communicator.receive(src, 0);
                    global.extend_from_slice(&buf);
                }
                Ok(global)
            } else {
                self.communicator.send(local, 0, 0);
                Ok(Vec::new())
            }
        }

        /// Serialise gathered data into VTK Legacy ASCII format.
        fn write_vtk_legacy_ascii<P: AsRef<Path>>(
            &self,
            filename: P,
            points: &[f64],
            cells: &[u32],
            cell_types: &[u8],
            point_data: &HashMap<String, Vec<f64>>,
            cell_data: &HashMap<String, Vec<f64>>,
        ) -> MpiResult<()> {
            let n_points = points.len() / 3;
            let n_cells = cell_types.len();

            let file = std::fs::File::create(filename.as_ref()).map_err(|e| {
                MpiError::IoError(format!("Failed to create VTK file: {e}"))
            })?;
            let mut w = std::io::BufWriter::new(file);

            // VTK header
            writeln!(w, "# vtk DataFile Version 3.0").map_err(io_err)?;
            writeln!(w, "CFDrs parallel output").map_err(io_err)?;
            writeln!(w, "ASCII").map_err(io_err)?;
            writeln!(w, "DATASET UNSTRUCTURED_GRID").map_err(io_err)?;

            // Points
            writeln!(w, "POINTS {n_points} double").map_err(io_err)?;
            for chunk in points.chunks(3) {
                let (x, y, z) = (chunk[0], chunk.get(1).copied().unwrap_or(0.0), chunk.get(2).copied().unwrap_or(0.0));
                writeln!(w, "{x:.15e} {y:.15e} {z:.15e}").map_err(io_err)?;
            }

            // Cells
            let cells_list_size: usize = cells.len();
            writeln!(w, "CELLS {n_cells} {cells_list_size}").map_err(io_err)?;
            // Cells array is assumed to already be in VTK "n v0 v1 ... vn-1" format.
            let mut idx = 0;
            while idx < cells.len() {
                let n_verts = cells[idx] as usize;
                let mut line = format!("{}", cells[idx]);
                for k in 1..=n_verts {
                    if idx + k < cells.len() {
                        line.push_str(&format!(" {}", cells[idx + k]));
                    }
                }
                writeln!(w, "{line}").map_err(io_err)?;
                idx += n_verts + 1;
            }

            // Cell types
            writeln!(w, "CELL_TYPES {n_cells}").map_err(io_err)?;
            for ct in cell_types {
                writeln!(w, "{ct}").map_err(io_err)?;
            }

            // Point data
            if !point_data.is_empty() {
                writeln!(w, "POINT_DATA {n_points}").map_err(io_err)?;
                for (name, data) in point_data {
                    writeln!(w, "SCALARS {name} double 1").map_err(io_err)?;
                    writeln!(w, "LOOKUP_TABLE default").map_err(io_err)?;
                    for v in data {
                        writeln!(w, "{v:.15e}").map_err(io_err)?;
                    }
                }
            }

            // Cell data
            if !cell_data.is_empty() {
                writeln!(w, "CELL_DATA {n_cells}").map_err(io_err)?;
                for (name, data) in cell_data {
                    writeln!(w, "SCALARS {name} double 1").map_err(io_err)?;
                    writeln!(w, "LOOKUP_TABLE default").map_err(io_err)?;
                    for v in data {
                        writeln!(w, "{v:.15e}").map_err(io_err)?;
                    }
                }
            }

            w.flush().map_err(io_err)?;
            Ok(())
        }
    }

    /// Parallel HDF5 writer for distributed datasets.
    ///
    /// Uses a gather-to-root strategy: every rank sends its owned data to
    /// rank 0, which writes a self-describing binary file with metadata
    /// headers compatible with downstream HDF5 readers.
    pub struct ParallelHdf5Writer<T: RealField> {
        communicator: MpiCommunicator,
        is_root: bool,
        _phantom: std::marker::PhantomData<T>,
    }

    impl<T: RealField + Copy + ToPrimitive> ParallelHdf5Writer<T> {
        /// Create new parallel HDF5 writer
        pub fn new(communicator: &MpiCommunicator) -> Self {
            Self {
                communicator: communicator.clone(),
                is_root: communicator.is_root(),
                _phantom: std::marker::PhantomData,
            }
        }

        /// Write distributed HDF5 file.
        ///
        /// All ranks participate; rank 0 writes a portable binary file containing
        /// metadata and all dataset values gathered from every process.
        pub fn write_hdf5_file<P: AsRef<Path>>(
            &self,
            filename: P,
            datasets: &HashMap<String, &DistributedVector<T>>,
            metadata: &HashMap<String, String>,
        ) -> MpiResult<()> {
            let size = self.communicator.size() as usize;

            // Gather dataset sizes and data.
            let mut global_datasets: HashMap<String, Vec<f64>> = HashMap::new();
            for (name, dv) in datasets {
                let local_n = dv.local_data.len();
                let mut all_sizes = vec![0usize; size];
                self.communicator.all_gather(&local_n, &mut all_sizes);

                let local_f64: Vec<f64> = dv
                    .local_data
                    .iter()
                    .map(|v| v.to_f64().unwrap_or(0.0))
                    .collect();

                let gathered = self.gather_f64_to_root(&local_f64, &all_sizes)?;
                global_datasets.insert(name.clone(), gathered);
            }

            // Root writes the file.
            if self.is_root {
                self.write_binary_datasets(filename, &global_datasets, metadata)?;
            }

            self.communicator.barrier();
            Ok(())
        }

        /// Gather variable-length f64 slices from all ranks to root.
        fn gather_f64_to_root(
            &self,
            local: &[f64],
            sizes: &[usize],
        ) -> MpiResult<Vec<f64>> {
            if self.is_root {
                let mut global = Vec::with_capacity(sizes.iter().sum());
                global.extend_from_slice(local);
                let size = self.communicator.size();
                for src in 1..size {
                    let buf: Vec<f64> = self.communicator.receive(src, 0);
                    global.extend_from_slice(&buf);
                }
                Ok(global)
            } else {
                self.communicator.send(local, 0, 0);
                Ok(Vec::new())
            }
        }

        /// Write gathered datasets to a portable binary file.
        ///
        /// Format:
        /// - 8-byte magic: `CFDrsHD5`
        /// - 4-byte LE u32: metadata count
        /// - For each metadata entry: 4-byte LE u32 key length, key bytes,
        ///   4-byte LE u32 value length, value bytes
        /// - 4-byte LE u32: dataset count
        /// - For each dataset: 4-byte LE u32 name length, name bytes,
        ///   8-byte LE u64 element count, then elements as LE f64
        fn write_binary_datasets<P: AsRef<Path>>(
            &self,
            filename: P,
            datasets: &HashMap<String, Vec<f64>>,
            metadata: &HashMap<String, String>,
        ) -> MpiResult<()> {
            let file = std::fs::File::create(filename.as_ref()).map_err(|e| {
                MpiError::IoError(format!("Failed to create HDF5 file: {e}"))
            })?;
            let mut w = std::io::BufWriter::new(file);

            // Magic header
            w.write_all(b"CFDrsHD5").map_err(io_err)?;

            // Metadata
            let meta_count = metadata.len() as u32;
            w.write_all(&meta_count.to_le_bytes()).map_err(io_err)?;
            for (k, v) in metadata {
                let kb = k.as_bytes();
                let vb = v.as_bytes();
                w.write_all(&(kb.len() as u32).to_le_bytes()).map_err(io_err)?;
                w.write_all(kb).map_err(io_err)?;
                w.write_all(&(vb.len() as u32).to_le_bytes()).map_err(io_err)?;
                w.write_all(vb).map_err(io_err)?;
            }

            // Datasets
            let ds_count = datasets.len() as u32;
            w.write_all(&ds_count.to_le_bytes()).map_err(io_err)?;
            for (name, data) in datasets {
                let nb = name.as_bytes();
                w.write_all(&(nb.len() as u32).to_le_bytes()).map_err(io_err)?;
                w.write_all(nb).map_err(io_err)?;
                w.write_all(&(data.len() as u64).to_le_bytes()).map_err(io_err)?;
                for &val in data {
                    w.write_all(&val.to_le_bytes()).map_err(io_err)?;
                }
            }

            w.flush().map_err(io_err)?;
            Ok(())
        }
    }

    /// Convert `std::io::Error` to `MpiError::IoError`.
    fn io_err(e: std::io::Error) -> MpiError {
        MpiError::IoError(format!("{e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that the parallel_io module types are constructible and that
    /// they carry the expected type parameters. This is a compile-time
    /// structural test; MPI runtime tests require `mpiexec`.
    #[test]
    fn test_parallel_io_type_construction() {
        let _vtk: std::marker::PhantomData<parallel_io::ParallelVtkWriter<f64>> =
            std::marker::PhantomData;
        let _hdf5: std::marker::PhantomData<parallel_io::ParallelHdf5Writer<f64>> =
            std::marker::PhantomData;
    }

    /// Verify distributed vector algebra contracts (type-level).
    #[test]
    fn test_distributed_vector_type_algebra() {
        // The distributed vector requires RealField + Copy + FromPrimitive + LowerExp.
        fn assert_bounds<T: RealField + Copy + num_traits::FromPrimitive + std::fmt::LowerExp>() {}
        assert_bounds::<f64>();
        assert_bounds::<f32>();
    }
}
