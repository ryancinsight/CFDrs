use super::traits::{DistributedLinearOperator, Preconditioner};
use super::vector::DistributedVector;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::decomposition::LocalSubdomain;
use nalgebra::RealField;
use num_traits::FromPrimitive;

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
            nx_local: 0,
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
    ) -> crate::compute::mpi::error::MpiResult<DistributedVector<T>> {
        let mut x = x0.clone();
        let mut residual_norm;

        // Initial residual: r = b - A*x0
        self.operator.apply(&x, &mut self.residual)?;
        self.residual.axpy(-T::one(), b);
        self.residual.scale(-T::one()); // r = b - A*x0

        residual_norm = self.residual.norm()?;
        let mut iter = 0;

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
    ) -> crate::compute::mpi::error::MpiResult<()> {
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
