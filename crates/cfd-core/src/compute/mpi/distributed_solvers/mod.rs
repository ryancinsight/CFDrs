//! Distributed linear solvers for MPI-parallel CFD computations.
//!
//! This module provides scalable iterative solvers that operate on
//! domain-decomposed meshes using MPI for inter-process communication. The
//! solver stack follows a preconditioned Krylov framework:
//!
//! 1. **`DistributedVector`** — MPI-aware vector with global reductions
//! 2. **`DistributedLaplacian2D`** — matrix-free discrete Laplacian operator
//! 3. **`BlockJacobiPreconditioner`** + **`AdditiveSchwarzPreconditioner`** — domain-decomposition
//!    preconditioners
//! 4. **`DistributedGMRES`** — restarted GMRES with left preconditioning
//! 5. **`parallel_io`** — gather-to-root VTK and HDF5 writers
//!
//! # Architecture
//!
//! ```text
//!  Vector ←─ Traits (Operator, Preconditioner)
//!             ├──► Laplacian2D   (impl Operator)
//!             ├──► BlockJacobi   (impl Preconditioner)
//!             ├──► AdditiveSchwarz (impl Preconditioner)
//!             └──► GMRES          (uses Operator + Preconditioner)
//!  ParallelIO ←── Vector
//! ```
//!
//! # Example
//!
//! ```ignore
//! // Create distributed operator and preconditioner
//! let laplacian = DistributedLaplacian2D::new(&decomp, &comm, dx, dy)?;
//! let precond = BlockJacobiPreconditioner::new(&laplacian, &decomp, &comm)?;
//!
//! // Create GMRES solver with 30-dim Krylov subspace
//! let mut gmres = DistributedGMRES::new(laplacian, precond, &comm, 30);
//!
//! // Solve: A x = b
//! let solution = gmres.solve(&b, &x0, 1e-10, 1000)?;
//! ```

/// Distributed linear operator and preconditioner traits.
pub mod traits;

/// MPI-aware distributed vector with global reductions.
pub mod vector;

/// Matrix-free 2D Laplacian operator for domain-decomposed grids.
pub mod laplacian;

/// Block Jacobi and additive Schwarz preconditioners.
pub mod preconditioners;

/// Distributed GMRES iterative solver.
pub mod gmres;

/// Parallel I/O writers (VTK, HDF5) using gather-to-root strategy.
pub mod parallel_io;

pub use gmres::DistributedGMRES;
pub use laplacian::DistributedLaplacian2D;
pub use parallel_io::{ParallelHdf5Writer, ParallelVtkWriter};
pub use preconditioners::{AdditiveSchwarzPreconditioner, BlockJacobiPreconditioner};
pub use traits::{DistributedLinearOperator, Preconditioner};
pub use vector::DistributedVector;

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
        fn assert_bounds<T: nalgebra::RealField + Copy + num_traits::FromPrimitive + std::fmt::LowerExp>(
        ) {
        }
        assert_bounds::<f64>();
        assert_bounds::<f32>();
    }
}
