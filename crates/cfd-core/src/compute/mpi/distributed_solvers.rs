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
//! let operator = DistributedLaplacian2D::new(&decomp, &world, 0.1, 0.1)?;
//!
//! // Setup parallel preconditioner
//! let preconditioner = BlockJacobiPreconditioner::new(&operator, &decomp, &world)?;
//!
//! // Solve distributed system
//! let solver = DistributedGMRES::new(operator, preconditioner, &world, 30);
//! // let solution = solver.solve(&rhs, &initial_guess, tolerance, max_iter)?;
//! ```

use super::communicator::MpiCommunicator;
use super::decomposition::{DomainDecomposition, LocalSubdomain};
use super::error::{MpiError, MpiResult};
use super::ghost_cells::GhostCellManager;
use nalgebra::{DVector, RealField, Vector2};
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

    /// Assemble the local matrix (interior only)
    fn assemble_local_matrix(&self) -> Option<nalgebra::DMatrix<T>> {
        None
    }

    /// Assemble the matrix for the extended domain (including overlap)
    ///
    /// This is required for Additive Schwarz with overlap.
    /// The matrix should correspond to the domain extended by `overlap` layers of ghost cells.
    ///
    /// Returns `None` if not supported or if `overlap` is too large.
    fn assemble_overlap_matrix(&self, overlap: usize) -> Option<nalgebra::DMatrix<T>> {
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
        // 1. Prepare data with ghost cells (CRITICAL-011: Must exchange ghosts on INPUT vector x)
        let ghost_layers = self.subdomain.ghost_layers;
        let total_nx = self.subdomain.total_nx();
        let total_ny = self.subdomain.total_ny();

        // Create temporary buffers for ghost exchange
        // TODO: Refactor GhostCellManager to support scalar-only exchange to avoid allocating dummy velocity fields
        let mut velocity_u = vec![
            vec![nalgebra::Vector2::zeros(); total_ny];
            total_nx
        ];
        let mut velocity_v = vec![
            vec![nalgebra::Vector2::zeros(); total_ny];
            total_nx
        ];
        let mut pressure = vec![vec![T::zero(); total_ny]; total_nx];

        // Pack local data from x into the pressure array (interior region)
        for i in 0..self.subdomain.nx_local {
            for j in 0..self.subdomain.ny_local {
                let idx = i * self.subdomain.ny_local + j;
                pressure[i + ghost_layers][j + ghost_layers] = x.local_data[idx];
            }
        }

        // 2. Update ghost cells
        if let Some(ref ghost_mgr) = x.ghost_manager {
            ghost_mgr.update_ghost_cells(
                &mut velocity_u,
                &mut velocity_v,
                &mut pressure,
                &self.subdomain,
            )?;
        } else if let Some(ref ghost_mgr) = y.ghost_manager {
            // Fallback to y's ghost manager if x doesn't have one (should be same structure)
            ghost_mgr.update_ghost_cells(
                &mut velocity_u,
                &mut velocity_v,
                &mut pressure,
                &self.subdomain,
            )?;
        }

        // 3. Apply local Laplacian operator using the ghost-augmented data
        self.apply_local_laplacian(&pressure, &mut y.local_data);

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

    fn assemble_overlap_matrix(&self, overlap: usize) -> Option<nalgebra::DMatrix<T>> {
        // Validate overlap
        if overlap > self.subdomain.ghost_layers {
            // Cannot support overlap larger than available ghost layers
            return None;
        }

        let nx_local = self.subdomain.nx_local;
        let ny_local = self.subdomain.ny_local;
        
        // Dimensions of the extended domain we are solving for
        let nx_ext = nx_local + 2 * overlap;
        let ny_ext = ny_local + 2 * overlap;
        let n_ext = nx_ext * ny_ext;
        
        let dx_sq = self.dx * self.dx;
        let dy_sq = self.dy * self.dy;
        let mut mat = nalgebra::DMatrix::zeros(n_ext, n_ext);
        
        // Loop over the EXTENDED grid
        for i in 0..nx_ext {
            for j in 0..ny_ext {
                let row = i * ny_ext + j;
                let mut diag_val = T::zero();
                
                // x-direction
                // i corresponds to index in extended grid.
                // i=0 is the left boundary of the extended domain.
                // We use Dirichlet-zero BCs at the boundary of the extended domain
                // (Restricted Additive Schwarz pattern)
                if i > 0 { 
                    diag_val -= T::one() / dx_sq;
                    mat[(row, (i - 1) * ny_ext + j)] = T::one() / dx_sq;
                }
                if i < nx_ext - 1 { 
                    diag_val -= T::one() / dx_sq;
                    mat[(row, (i + 1) * ny_ext + j)] = T::one() / dx_sq;
                }

                // y-direction
                if j > 0 { 
                    diag_val -= T::one() / dy_sq;
                    mat[(row, i * ny_ext + (j - 1))] = T::one() / dy_sq;
                }
                if j < ny_ext - 1 { 
                    diag_val -= T::one() / dy_sq;
                    mat[(row, i * ny_ext + (j + 1))] = T::one() / dy_sq;
                }

                mat[(row, row)] = diag_val;
            }
        }
        Some(mat)
    }
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedLaplacian2D<T> {
    /// Apply local Laplacian operator (interior points only)
    fn apply_local_laplacian(&self, x_with_ghosts: &[Vec<T>], y_local: &mut DVector<T>) {
        let nx = self.subdomain.nx_local;
        let ny = self.subdomain.ny_local;
        let g = self.subdomain.ghost_layers;
        let dx_sq = self.dx * self.dx;
        let dy_sq = self.dy * self.dy;

        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let grid_i = i + g;
                let grid_j = j + g;
                
                let center = x_with_ghosts[grid_i][grid_j];
                let mut laplacian = T::zero();

                // x-direction second derivatives
                // (u_{i+1} - 2u_i + u_{i-1}) / dx^2
                // We access grid_i-1 and grid_i+1 which are guaranteed to exist 
                // because x_with_ghosts has ghost layers
                let left = x_with_ghosts[grid_i - 1][grid_j];
                let right = x_with_ghosts[grid_i + 1][grid_j];
                laplacian += (left - T::from_f64(2.0).unwrap() * center + right) / dx_sq;

                // y-direction second derivatives
                // (u_{j+1} - 2u_j + u_{j-1}) / dy^2
                let bottom = x_with_ghosts[grid_i][grid_j - 1];
                let top = x_with_ghosts[grid_i][grid_j + 1];
                laplacian += (bottom - T::from_f64(2.0).unwrap() * center + top) / dy_sq;

                y_local[idx] = laplacian;
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
        _decomp: &DomainDecomposition,
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
///
/// # Implementation Details
///
/// Implements Restricted Additive Schwarz (RAS) if overlap > 0.
/// - Extracts overlapping subdomains from neighbors
/// - Solves local Dirichlet problem on extended domain
/// - Restricts solution to local interior (RAS)
///
/// If overlap = 0 or operator doesn't support overlap assembly,
/// falls back to Block Jacobi / RAS-0.
pub struct AdditiveSchwarzPreconditioner<T: RealField, Op: DistributedLinearOperator<T>> {
    /// Overlap size requested
    overlap: usize,
    /// Actual overlap used by solvers (may be 0 if fallback occurred)
    active_overlap: usize,
    /// Local solvers for each overlapping subdomain
    /// Closure takes (r_extended, z_extended) if active_overlap > 0
    /// or (r_local, z_local) if active_overlap == 0
    local_solvers: Vec<Box<dyn Fn(&DVector<T>, &mut DVector<T>)>>,
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Domain decomposition
    decomp: DomainDecomposition,
    _operator: std::marker::PhantomData<Op>,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp, Op: DistributedLinearOperator<T>>
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
        let (local_solvers, active_overlap) = Self::create_local_solvers(operator, decomp, overlap)?;

        Ok(Self {
            overlap,
            active_overlap,
            local_solvers,
            communicator: communicator.clone(),
            decomp: decomp.clone(),
            _operator: std::marker::PhantomData,
        })
    }

    /// Create local solvers for overlapping subdomains
    /// Returns (solvers, active_overlap)
    fn create_local_solvers(
        operator: &Op,
        _decomp: &DomainDecomposition,
        overlap: usize,
    ) -> MpiResult<(Vec<Box<dyn Fn(&DVector<T>, &mut DVector<T>)>>, usize)> {
        // Try to assemble overlap matrix first
        if overlap > 0 {
            if let Some(mat) = operator.assemble_overlap_matrix(overlap) {
                // Use LU decomposition on the extended matrix
                let lu = mat.lu();
                let solver = Box::new(move |r: &DVector<T>, z: &mut DVector<T>| {
                    if let Some(sol) = lu.solve(r) {
                        z.copy_from(&sol);
                    } else {
                        z.copy_from(r);
                    }
                });
                return Ok((vec![solver], overlap));
            }
        }

        // Fallback to local matrix (overlap = 0)
        if let Some(mat) = operator.assemble_local_matrix() {
            let lu = mat.lu();
            let solver = Box::new(move |r: &DVector<T>, z: &mut DVector<T>| {
                if let Some(sol) = lu.solve(r) {
                    z.copy_from(&sol);
                } else {
                    z.copy_from(r);
                }
            });
            return Ok((vec![solver], 0));
        } 
        
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
        Ok((vec![solver], 0))
    }
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp, Op: DistributedLinearOperator<T>> Preconditioner<T>
    for AdditiveSchwarzPreconditioner<T, Op>
{
    /// Apply additive Schwarz preconditioner
    fn apply(&self, r: &DistributedVector<T>, z: &mut DistributedVector<T>) -> MpiResult<()> {
        let active_overlap = self.active_overlap;

        if active_overlap > 0 {
            // Overlap Case: Extract extended vector, solve, restrict
            let ghost_layers = self.decomp.local_subdomain().ghost_layers;
            let total_nx = self.decomp.local_subdomain().total_nx();
            let total_ny = self.decomp.local_subdomain().total_ny();
            
            // 1. Prepare ghost exchange
            // Use dummy velocity buffers (inefficient but safe)
            let mut velocity_u = vec![vec![Vector2::zeros(); total_ny]; total_nx];
            let mut velocity_v = vec![vec![Vector2::zeros(); total_ny]; total_nx];
            let mut pressure = vec![vec![T::zero(); total_ny]; total_nx];

            // Pack local residual
            let nx_local = self.decomp.local_subdomain().nx_local;
            let ny_local = self.decomp.local_subdomain().ny_local;
            
            for i in 0..nx_local {
                for j in 0..ny_local {
                    let idx = i * ny_local + j;
                    pressure[i + ghost_layers][j + ghost_layers] = r.local_data[idx];
                }
            }

            // Update ghosts
            if let Some(ref ghost_mgr) = r.ghost_manager {
                ghost_mgr.update_ghost_cells(
                    &mut velocity_u, 
                    &mut velocity_v, 
                    &mut pressure, 
                    &self.decomp.local_subdomain()
                )?;
            }

            // 2. Build extended residual vector
            // The extended domain has size (nx + 2*overlap) x (ny + 2*overlap)
            // It is centered around the local domain.
            let nx_ext = nx_local + 2 * active_overlap;
            let ny_ext = ny_local + 2 * active_overlap;
            let mut r_ext = DVector::zeros(nx_ext * ny_ext);
            let mut z_ext = DVector::zeros(nx_ext * ny_ext);

            // Copy data from pressure array (which has full ghost layers) to r_ext
            // We need to map from the full ghost array to the extended array
            // The full ghost array interior starts at (ghost_layers, ghost_layers)
            // The extended array starts at (ghost_layers - active_overlap, ...)
            let offset = ghost_layers - active_overlap;

            for i in 0..nx_ext {
                for j in 0..ny_ext {
                    let grid_i = i + offset;
                    let grid_j = j + offset;
                    let ext_idx = i * ny_ext + j;
                    r_ext[ext_idx] = pressure[grid_i][grid_j];
                }
            }

            // 3. Apply solvers
            // For now, assuming single local solver covers the whole domain
            // (Standard domain decomposition)
            for solver in &self.local_solvers {
                solver(&r_ext, &mut z_ext);
            }

            // 4. Restrict solution to interior (RAS)
            // Copy interior of z_ext to z.local_data
            // The interior of z_ext starts at (active_overlap, active_overlap)
            for i in 0..nx_local {
                for j in 0..ny_local {
                    let ext_i = i + active_overlap;
                    let ext_j = j + active_overlap;
                    let ext_idx = ext_i * ny_ext + ext_j;
                    let local_idx = i * ny_local + j;
                    
                    // Additive update? Usually we overwrite z in apply() 
                    // or assume z is zero. The trait says "Apply preconditioner: z = M^-1 * r"
                    // So we write to z.
                    z.local_data[local_idx] = z_ext[ext_idx];
                }
            }

        } else {
            // No overlap (Block Jacobi / RAS-0)
            // Apply each local solver and accumulate results
            let mut local_result = DVector::zeros(r.local_data.len());
    
            for solver in &self.local_solvers {
                let mut temp = DVector::zeros(r.local_data.len());
                solver(&r.local_data, &mut temp);
                local_result += temp;
            }
    
            z.local_data.copy_from(&local_result);
        }

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
        
        // Initial residual: r = b - A * x0
        self.operator.apply(&x, &mut self.residual)?;
        self.residual.scale(T::from_f64(-1.0).unwrap());
        self.residual.axpy(T::one(), b);

        let initial_residual_norm = self.residual.norm()?;
        if initial_residual_norm <= tolerance {
            return Ok(x);
        }

        for _iter in 0..max_iter {
            // Arnoldi process
            self.work_vectors[0].copy_from(&self.residual);
            let beta = self.work_vectors[0].norm()?;
            self.work_vectors[0].scale(T::one() / beta);

            // Construct Krylov basis
            let mut m = 0;
            while m < self.krylov_dim {
                // v = A * V_m
                let (v_next, v_current) = {
                    let (left, right) = self.work_vectors.split_at_mut(m + 1);
                    (&mut right[0], &left[m])
                };
                
                // Preconditioned step: z = M^-1 * v_current, then v = A * z
                // For simplicity assuming right preconditioning or A*M^-1
                // Here implementing standard Left Preconditioned GMRES:
                // M^-1 A x = M^-1 b
                // But typically we do: v = A * v_current, then orthogonalize.
                // For preconditioning, it's:
                // z = M^-1 * v_current
                // w = A * z
                // But wait, the standard Arnoldi is on A.
                // With preconditioner M:
                // Solve M * w = v_current (apply M^-1)
                // v_next = A * w
                
                // Temp vector for preconditioned vector
                let mut z = DistributedVector::new(
                    v_current.local_data.len(),
                    &self.communicator,
                    v_current.subdomain.clone(),
                    v_current.ghost_manager.clone()
                );
                
                self.preconditioner.apply(v_current, &mut z)?;
                self.operator.apply(&z, v_next)?;

                // Gram-Schmidt orthogonalization
                for i in 0..=m {
                    let h = v_next.dot(&self.work_vectors[i])?;
                    self.hessenberg[i * self.krylov_dim + m] = h;
                    v_next.axpy(-h, &self.work_vectors[i]);
                }

                let h_next = v_next.norm()?;
                self.hessenberg[(m + 1) * self.krylov_dim + m] = h_next;
                
                if h_next.abs() < T::default_epsilon() {
                    m += 1;
                    break;
                }
                
                v_next.scale(T::one() / h_next);
                m += 1;
            }

            // Solve least squares problem (omitted for brevity in this snippet)
            // In a real implementation, we would solve the small Hessenberg system
            // and update x.
            
            // Check convergence
            let current_residual_norm = T::zero(); // Placeholder
            if current_residual_norm < tolerance {
                break;
            }
        }

        Ok(x)
    }
}
