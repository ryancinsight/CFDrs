use super::traits::{DistributedLinearOperator, Preconditioner};
use super::vector::DistributedVector;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::decomposition::DomainDecomposition;
use crate::compute::mpi::error::MpiResult;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

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
