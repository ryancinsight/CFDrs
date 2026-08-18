//! Block Preconditioners for Saddle-Point Systems
//!
//! Provides specialized preconditioners for incompressible flow systems:
//! ```text
//! [ A   B^T ] [ u ]   [ f ]
//! [ B   0   ] [ p ] = [ g ]
//! ```
//!
//! where:
//! - A: Momentum matrix (viscous + convection)
//! - B: Divergence operator
//! - B^T: Gradient operator
//!
//! # References
//!
//! 1. Elman, Silvester & Wathen (2005): "Finite Elements and Fast Iterative Solvers"
//! 2. Benzi, Golub & Liesen (2005): "Numerical solution of saddle point problems"
//! 3. Murphy, Golub & Wathen (2000): "A Note on Preconditioning for Indefinite Linear Systems"

use crate::linear_solver::Preconditioner;
use crate::sparse::SparseMatrix;
use cfd_core::error::Result;
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto::LetoError;
use leto_ops::{
    factor_symbolic, CscMatrix, OwnedNumericLu, RealScalar as LetoRealScalar, Scalar as LetoScalar,
    SparseLuSolver, SymbolicLu,
};

use crate::sparse::SparseMatrixBuilder;

#[inline]
fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: usize fits in u64");
    <T as FloatElement>::from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

#[inline]
fn diagonal_epsilon<T: FloatElement>() -> T {
    <T as FloatElement>::from_f64(1e-14)
}

#[inline]
fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

fn validate_vector_len<T>(
    name: &str,
    vector: &Array1<T>,
    expected: usize,
) -> std::result::Result<(), LetoError> {
    let actual = vector_len(vector);
    if actual != expected {
        return Err(LetoError::InvalidInput(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

/// Add provider-pivot-scale diagonal entries to structurally empty or
/// numerically zero preconditioner rows.
///
/// Component extraction can leave an isolated velocity or pressure row even
/// when the full saddle system is valid. A direct factorization must reject
/// that singular block; the block preconditioner instead represents the row by
/// the provider's pivot-scale identity, which preserves a bounded solve while
/// leaving the assembled operator unchanged.
fn stabilize_preconditioner_diagonal<T>(
    matrix: &SparseMatrix<T>,
    pivot_tolerance: f64,
) -> Result<SparseMatrix<T>>
where
    T: RealField + FloatElement + Copy + LetoScalar,
{
    let n = matrix.nrows();
    let mut builder = SparseMatrixBuilder::new(n, matrix.ncols());
    let pivot_scale = <T as FloatElement>::from_f64(pivot_tolerance);

    for row_index in 0..n {
        let row = matrix.row(row_index);
        let mut diagonal = <T as NumericElement>::ZERO;
        let mut row_scale = <T as NumericElement>::ZERO;
        for (&column, &value) in row.col_indices().iter().zip(row.values()) {
            builder.add_entry(row_index, column, value)?;
            let magnitude = NumericElement::abs(value);
            if magnitude > row_scale {
                row_scale = magnitude;
            }
            if column == row_index {
                diagonal += value;
            }
        }

        if NumericElement::abs(diagonal) <= diagonal_epsilon() {
            let mut regularization = row_scale * pivot_scale;
            if regularization < diagonal_epsilon() {
                regularization = diagonal_epsilon();
            }
            builder.add_entry(row_index, row_index, regularization)?;
        }
    }

    builder.build()
}

/// Build a pressure preconditioner from signed row-sum lumping of a Schur
/// approximation. Preserving the sign is required for the indefinite saddle
/// system; if row entries cancel, the raw diagonal supplies the signed scale.
fn lumped_pressure_preconditioner<T>(matrix: &SparseMatrix<T>) -> DiagonalPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoScalar,
{
    let n = matrix.nrows();
    let mut diag_inv = Array1::zeros([n]);
    for row_index in 0..n {
        let row_sum = matrix
            .row(row_index)
            .values()
            .iter()
            .copied()
            .fold(<T as NumericElement>::ZERO, |sum, value| sum + value);
        let scale = if NumericElement::abs(row_sum) > diagonal_epsilon() {
            row_sum
        } else {
            get_diagonal(matrix, row_index)
        };
        diag_inv[row_index] = if NumericElement::abs(scale) > diagonal_epsilon() {
            <T as NumericElement>::ONE / scale
        } else {
            <T as NumericElement>::ONE
        };
    }
    DiagonalPreconditioner { diag_inv }
}

/// Extract diagonal element from CSR matrix
fn get_diagonal<T: RealField + Copy + LetoScalar>(matrix: &SparseMatrix<T>, row: usize) -> T {
    let row_range = matrix.row(row);
    for (col_idx, &value) in row_range.col_indices().iter().zip(row_range.values()) {
        if *col_idx == row {
            return value;
        }
    }
    <T as NumericElement>::ZERO // Return zero if diagonal entry not found
}

/// Block diagonal preconditioner for saddle-point systems
///
/// P = [ A_inv    0   ]
///     [ 0      S_inv ]
///
/// where S ≈ B A^{-1} B^T (Schur complement approximation)
///
/// # Algorithm
///
/// For solving Px = b:
/// 1. Solve A u_tilde = f  (momentum block)
/// 2. Solve S p = g - B u_tilde  (pressure Schur complement)
/// 3. Solve A u = f - B^T p  (pressure correction)
///
/// # Theorem (Block-Diagonal Preconditioning of Saddle-Point Systems)
///
/// For the saddle-point system with SPD momentum block $A$ and full-rank $B$,
/// the exact block-diagonal preconditioner
/// $P = \mathrm{diag}(A,\;S)$ with $S = B A^{-1} B^T$ yields a preconditioned
/// system whose eigenvalues lie in $\{1\} \cup \bigl[\frac{1-\sqrt{5}}{2},\;\frac{1+\sqrt{5}}{2}\bigr]$,
/// guaranteeing GMRES convergence in at most 3 iterations.
///
/// When the Schur complement is approximated by $\tilde{S} \approx B\,\mathrm{diag}(A)^{-1}\,B^T$,
/// the spectral bounds degrade gracefully with the quality of the diagonal
/// approximation.
///
/// **Proof sketch**: The preconditioned matrix $P^{-1}\mathcal{A}$ satisfies
/// $(P^{-1}\mathcal{A})^3 - (P^{-1}\mathcal{A})^2 = 0$ when $S$ is exact,
/// giving eigenvalues from the characteristic polynomial $\lambda^2 - \lambda - 1 = 0$
/// plus $\lambda = 1$.
///
/// **Reference**: Murphy, Golub & Wathen (2000), Theorem 2.1;
/// Elman, Silvester & Wathen (2005), §6.2.
///
/// # Performance
///
/// - Complexity: O(nnz) per iteration (sparse matrix operations)
/// - Convergence: Significantly better than ILU for saddle-point systems
/// - Memory: Requires storing momentum matrix A and mass matrix M (for S approximation)
pub struct BlockDiagonalPreconditioner<T: RealField + FloatElement> {
    /// Momentum matrix approximation (e.g., diagonal or ILU of A)
    momentum_preconditioner: DiagonalPreconditioner<T>,
    /// Pressure Schur complement approximation (typically pressure mass matrix)
    pressure_preconditioner: DiagonalPreconditioner<T>,
    /// Velocity DOF count
    n_velocity: usize,
    /// Pressure DOF count
    n_pressure: usize,
}

/// Simple diagonal preconditioner (Jacobi)
pub struct DiagonalPreconditioner<T: RealField> {
    /// Inverse of diagonal entries
    diag_inv: Array1<T>,
}

impl<T: RealField + FloatElement + LetoScalar> DiagonalPreconditioner<T> {
    /// Create diagonal preconditioner from matrix diagonal
    pub fn new(matrix: &SparseMatrix<T>) -> Self {
        let n = matrix.nrows();
        let mut diag_inv = Array1::zeros([n]);

        for i in 0..n {
            let d = get_diagonal(matrix, i);
            if NumericElement::abs(d) > diagonal_epsilon() {
                diag_inv[i] = <T as NumericElement>::ONE / d;
            } else {
                // Singular diagonal entry, use identity
                diag_inv[i] = <T as NumericElement>::ONE;
            }
        }

        Self { diag_inv }
    }

    /// Apply preconditioner: x = M^{-1} b
    pub fn apply(&self, b: &Array1<T>) -> Result<Array1<T>> {
        validate_vector_len(
            "diagonal preconditioner input",
            b,
            vector_len(&self.diag_inv),
        )?;

        let mut x = Array1::zeros([vector_len(b)]);
        for idx in 0..vector_len(b) {
            x[idx] = b[idx] * self.diag_inv[idx];
        }
        Ok(x)
    }
}

impl<T: RealField + FloatElement + Copy + LetoScalar> BlockDiagonalPreconditioner<T> {
    /// Create block diagonal preconditioner
    ///
    /// # Arguments
    ///
    /// * `matrix` - Full saddle-point system matrix
    /// * `n_velocity` - Number of velocity DOFs
    /// * `n_pressure` - Number of pressure DOFs
    ///
    /// # Algorithm
    ///
    /// 1. Extract momentum block A (top-left n_velocity × n_velocity)
    /// 2. Extract gradient block B^T (top-right n_velocity × n_pressure)
    /// 3. Approximate Schur complement: S ≈ ν M_p (viscosity-scaled mass matrix)
    ///    - Better than diagonal: accounts for element coupling
    ///    - M_p_ij = ∫ φ_i φ_j dΩ (mass matrix assembly)
    ///    - Use: S^{-1} ≈ (1/ν) M_p^{-1}
    pub fn new(matrix: &SparseMatrix<T>, n_velocity: usize, n_pressure: usize) -> Result<Self> {
        if matrix.nrows() != n_velocity + n_pressure {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Matrix size mismatch: {} != {} + {}",
                matrix.nrows(),
                n_velocity,
                n_pressure
            )));
        }

        // Extract momentum block diagonal (A block)
        let mut momentum_diag = Array1::zeros([n_velocity]);
        for i in 0..n_velocity {
            let d = get_diagonal(matrix, i);
            momentum_diag[i] = if NumericElement::abs(d) > diagonal_epsilon() {
                d
            } else {
                <T as NumericElement>::ONE // Avoid division by zero
            };
        }

        // IMPROVED: Approximate Schur complement using viscosity-scaled mass matrix
        // For incompressible flow: S ≈ μ M_p where M_p is pressure mass matrix
        // Extract pressure block diagonal and off-diagonal contributions
        let mut pressure_diag = Array1::zeros([n_pressure]);

        for i in 0..n_pressure {
            let row_idx = n_velocity + i;

            // Sum row absolute values to approximate mass matrix diagonal
            // M_p_ii ≈ Σ_j |M_p_ij|  (row-sum lumping)
            let row = matrix.row(row_idx);
            let mut row_sum = <T as NumericElement>::ZERO;

            for (col_idx, &value) in row.col_indices().iter().zip(row.values()) {
                if *col_idx >= n_velocity {
                    // Only pressure-pressure block
                    row_sum += NumericElement::abs(value);
                }
            }

            if row_sum > diagonal_epsilon() {
                pressure_diag[i] = row_sum;
            } else {
                // Fallback: use viscosity-based scaling
                // Estimate from momentum block diagonal (average viscosity)
                let mut sum = <T as NumericElement>::ZERO;
                let mut count = 0usize;
                for idx in 0..n_velocity {
                    let value = NumericElement::abs(momentum_diag[idx]);
                    if value > diagonal_epsilon() {
                        sum += value;
                        count += 1;
                    }
                }

                let avg_momentum_diag = if count == 0 {
                    <T as NumericElement>::ONE
                } else {
                    sum / from_usize(count)
                };

                pressure_diag[i] = avg_momentum_diag;
            }
        }

        let mut momentum_diag_inv = Array1::zeros([n_velocity]);
        for idx in 0..n_velocity {
            let d = momentum_diag[idx];
            momentum_diag_inv[idx] = if NumericElement::abs(d) > diagonal_epsilon() {
                <T as NumericElement>::ONE / d
            } else {
                <T as NumericElement>::ONE
            };
        }

        let mut pressure_diag_inv = Array1::zeros([n_pressure]);
        for idx in 0..n_pressure {
            let d = pressure_diag[idx];
            pressure_diag_inv[idx] = if NumericElement::abs(d) > diagonal_epsilon() {
                <T as NumericElement>::ONE / d
            } else {
                <T as NumericElement>::ONE
            };
        }

        let momentum_preconditioner = DiagonalPreconditioner {
            diag_inv: momentum_diag_inv,
        };

        let pressure_preconditioner = DiagonalPreconditioner {
            diag_inv: pressure_diag_inv,
        };

        Ok(Self {
            momentum_preconditioner,
            pressure_preconditioner,
            n_velocity,
            n_pressure,
        })
    }

    /// Apply block diagonal preconditioner: x = P^{-1} b
    ///
    /// # Algorithm
    ///
    /// Split b = [f, g]^T into velocity and pressure parts:
    /// 1. u = A_inv * f  (momentum preconditioning)
    /// 2. p = S_inv * g  (pressure preconditioning)
    /// 3. Return [u, p]^T
    pub fn apply(&self, b: &Array1<T>) -> Result<Array1<T>> {
        let n_total = self.n_velocity + self.n_pressure;
        validate_vector_len("block diagonal preconditioner input", b, n_total)?;

        let mut x = Array1::zeros([n_total]);

        // Apply momentum block preconditioner
        for idx in 0..self.n_velocity {
            x[idx] = b[idx] * self.momentum_preconditioner.diag_inv[idx];
        }

        // Apply pressure block preconditioner
        for idx in 0..self.n_pressure {
            x[self.n_velocity + idx] =
                b[self.n_velocity + idx] * self.pressure_preconditioner.diag_inv[idx];
        }

        Ok(x)
    }
}

/// SIMPLE preconditioner (Semi-Implicit Method for Pressure-Linked Equations)
///
/// More sophisticated than block diagonal, uses coupling between u and p.
///
/// # Algorithm
///
/// For solving $\begin{bmatrix} A & G \\ D & C \end{bmatrix} \begin{bmatrix} u \\ p \end{bmatrix} = \begin{bmatrix} f \\ g \end{bmatrix}$:
///
/// 1. Solve $A u^* = f$ (momentum prediction via diagonal approximation)
/// 2. Solve $S p = g - D u^*$ (pressure Schur equation)
/// 3. Correct $u = u^* - \text{diag}(A)^{-1} G p$
///
/// where $S \approx D\,\text{diag}(A)^{-1} G$ (Schur complement).
///
/// The diagonal Schur approximation is a block-structured heuristic. Its
/// effectiveness depends on the momentum block, coupling signs, stabilization,
/// and pressure nullspace; convergence is verified by the enclosing Krylov
/// solver rather than guaranteed by this preconditioner alone.
///
/// # References
///
/// Patankar, S. V. (1980): "Numerical Heat Transfer and Fluid Flow"
pub struct SimplePreconditioner<T: RealField + FloatElement> {
    /// Momentum block preconditioner (diag(A)^{-1})
    momentum_inv: DiagonalPreconditioner<T>,
    /// Inverse diagonal of the Schur complement approximation
    schur_diag_inv: Array1<T>,
    /// Rows of the divergence block stored as (velocity index, value) pairs
    /// per pressure row, used for the `D u*` product.
    divergence_rows: Vec<Vec<(usize, T)>>,
    /// Columns of the gradient block stored as (velocity index, value) pairs
    /// per pressure column, used for the `G p` correction.
    gradient_columns: Vec<Vec<(usize, T)>>,
    n_velocity: usize,
    n_pressure: usize,
}

impl<T: RealField + FloatElement + Copy + LetoScalar> SimplePreconditioner<T> {
    /// Create SIMPLE preconditioner.
    ///
    /// Extracts the $A$, $D$, and $G$ sub-blocks from the full saddle-point
    /// matrix and builds the diagonal Schur complement
    /// $\text{diag}(C - D\,\text{diag}(A)^{-1}G)$.
    pub fn new(matrix: &SparseMatrix<T>, n_velocity: usize, n_pressure: usize) -> Result<Self> {
        let eps = diagonal_epsilon();

        // Extract momentum diagonal and its inverse
        let mut momentum_diag = Array1::zeros([n_velocity]);
        for i in 0..n_velocity {
            momentum_diag[i] = get_diagonal(matrix, i);
        }

        let mut momentum_diag_inv = Array1::zeros([n_velocity]);
        for idx in 0..n_velocity {
            let d = momentum_diag[idx];
            momentum_diag_inv[idx] = if NumericElement::abs(d) > eps {
                <T as NumericElement>::ONE / d
            } else {
                <T as NumericElement>::ONE
            };
        }

        let momentum_inv = DiagonalPreconditioner {
            diag_inv: momentum_diag_inv,
        };

        // Extract the coupling blocks independently. The continuity row may
        // be scaled by the formulation, so `G` is not reconstructed from `D`.
        let mut divergence_rows = Vec::with_capacity(n_pressure);
        for i in 0..n_pressure {
            let global_row = n_velocity + i;
            let row = matrix.row(global_row);
            let entries: Vec<(usize, T)> = row
                .col_indices()
                .iter()
                .zip(row.values())
                .filter(|(&c, _)| c < n_velocity)
                .map(|(&c, &v)| (c, v))
                .collect();
            divergence_rows.push(entries);
        }
        let mut gradient_columns = vec![Vec::new(); n_pressure];
        for velocity_row in 0..n_velocity {
            let row = matrix.row(velocity_row);
            for (&column, &value) in row.col_indices().iter().zip(row.values()) {
                if let Some(pressure_column) = column
                    .checked_sub(n_velocity)
                    .filter(|&index| index < n_pressure)
                {
                    gradient_columns[pressure_column].push((velocity_row, value));
                }
            }
        }

        // Compute diag(C - D diag(A)^-1 G) for the actual assembled coupling
        // signs. This remains valid when the continuity row is normalized and
        // retains the pressure stabilization block C.
        let mut schur_diag_inv = Array1::zeros([n_pressure]);
        for i in 0..n_pressure {
            let mut s_ii = get_diagonal(matrix, n_velocity + i);
            for &(velocity_index, divergence_value) in &divergence_rows[i] {
                if let Some(&(_, gradient_value)) = gradient_columns[i]
                    .iter()
                    .find(|&&(index, _)| index == velocity_index)
                {
                    s_ii -=
                        divergence_value * momentum_inv.diag_inv[velocity_index] * gradient_value;
                }
            }
            schur_diag_inv[i] = if NumericElement::abs(s_ii) > eps {
                <T as NumericElement>::ONE / s_ii
            } else {
                <T as NumericElement>::ONE
            };
        }

        Ok(Self {
            momentum_inv,
            schur_diag_inv,
            divergence_rows,
            gradient_columns,
            n_velocity,
            n_pressure,
        })
    }

    /// Apply SIMPLE preconditioner with momentum-pressure coupling.
    ///
    /// Given $b = [f, g]^T$:
    /// 1. $u^* = \text{diag}(A)^{-1} f$
    /// 2. $p   = S^{-1}(g - D u^*)$
    /// 3. $u   = u^* - \text{diag}(A)^{-1} G p$
    pub fn apply(&self, b: &Array1<T>) -> Result<Array1<T>> {
        let n_total = self.n_velocity + self.n_pressure;
        validate_vector_len("SIMPLE preconditioner input", b, n_total)?;

        let mut x = Array1::zeros([n_total]);

        // Step 1: Momentum prediction u* = diag(A)^{-1} f
        let mut u_star = Array1::zeros([self.n_velocity]);
        for idx in 0..self.n_velocity {
            u_star[idx] = b[idx] * self.momentum_inv.diag_inv[idx];
        }

        // Step 2: Pressure correction p = S^{-1} (g - D u*)
        let mut rhs_p = Array1::zeros([self.n_pressure]);
        for idx in 0..self.n_pressure {
            rhs_p[idx] = b[self.n_velocity + idx];
        }
        for i in 0..self.n_pressure {
            let mut b_u = <T as NumericElement>::ZERO;
            for &(velocity_index, divergence_value) in &self.divergence_rows[i] {
                b_u += divergence_value * u_star[velocity_index];
            }
            rhs_p[i] -= b_u;
        }
        let mut p = Array1::zeros([self.n_pressure]);
        for idx in 0..self.n_pressure {
            p[idx] = rhs_p[idx] * self.schur_diag_inv[idx];
        }

        // Step 3: Velocity correction u = u* - diag(A)^{-1} G p
        let mut u_corrected = u_star;
        for i in 0..self.n_pressure {
            for &(velocity_index, gradient_value) in &self.gradient_columns[i] {
                u_corrected[velocity_index] -=
                    self.momentum_inv.diag_inv[velocity_index] * gradient_value * p[i];
            }
        }

        for idx in 0..self.n_velocity {
            x[idx] = u_corrected[idx];
        }
        for idx in 0..self.n_pressure {
            x[self.n_velocity + idx] = p[idx];
        }

        Ok(x)
    }
}

impl<T> Preconditioner<T> for BlockDiagonalPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoScalar,
{
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> std::result::Result<(), LetoError> {
        validate_vector_len(
            "block diagonal preconditioner input",
            r,
            self.n_velocity + self.n_pressure,
        )?;
        validate_vector_len("block diagonal preconditioner output", z, vector_len(r))?;
        for idx in 0..self.n_velocity {
            z[idx] = r[idx] * self.momentum_preconditioner.diag_inv[idx];
        }
        for idx in 0..self.n_pressure {
            let pressure_idx = self.n_velocity + idx;
            z[pressure_idx] = r[pressure_idx] * self.pressure_preconditioner.diag_inv[idx];
        }
        Ok(())
    }
}

impl<T> Preconditioner<T> for SimplePreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoScalar,
{
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> std::result::Result<(), LetoError> {
        validate_vector_len("SIMPLE preconditioner output", z, vector_len(r))?;
        let result = self
            .apply(r)
            .map_err(|e| LetoError::InvalidInput(e.to_string()))?;
        for idx in 0..vector_len(z) {
            z[idx] = result[idx];
        }
        Ok(())
    }
}

/// Number of velocity components in the three-dimensional Taylor–Hood system.
const VELOCITY_COMPONENTS: usize = 3;

/// Component-block sparse-LU preconditioner for a three-dimensional saddle system.
///
/// The velocity block in CFDrs is component-major. When the assembled operator
/// has no cross-component velocity entries, its momentum block is the direct
/// sum of three scalar CSR operators. This preconditioner factors each scalar
/// block once with the provider-owned Leto sparse LU implementation and combines
/// those solves with the independently assembled SIMPLE pressure correction.
///
/// The constructor rejects cross-component entries. Applying only diagonal
/// component blocks to a coupled operator would change the preconditioner
/// contract while appearing to succeed, so such systems fall through to the
/// general saddle-point tiers.
pub struct ComponentBlockPreconditioner<T: RealField + FloatElement + LetoRealScalar> {
    momentum_blocks: Vec<OwnedNumericLu<T>>,
    pressure_block: PressureBlock<T>,
    simple: SimplePreconditioner<T>,
    component_size: usize,
    n_velocity: usize,
    n_pressure: usize,
}

enum PressureBlock<T: RealField + FloatElement + LetoRealScalar> {
    LumpedDiagonal(DiagonalPreconditioner<T>),
}

/// Cached provider symbolic factors for an unchanged reduced FEM topology.
#[derive(Debug, Clone)]
pub(crate) struct ComponentBlockPattern {
    component_size: usize,
    n_velocity: usize,
    n_pressure: usize,
    row_ptr: Vec<usize>,
    col_indices: Vec<usize>,
    symbols: Vec<SymbolicLu>,
}

impl ComponentBlockPattern {
    fn matches<T: RealField + Copy + LetoScalar>(
        &self,
        matrix: &SparseMatrix<T>,
        n_velocity: usize,
        n_pressure: usize,
    ) -> bool {
        self.component_size * VELOCITY_COMPONENTS == n_velocity
            && self.n_velocity == n_velocity
            && self.n_pressure == n_pressure
            && self.row_ptr == matrix.row_ptr()
            && self.col_indices == matrix.col_indices()
    }
}

impl<T> ComponentBlockPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoRealScalar,
{
    /// Factor the independent velocity component blocks of `matrix`.
    ///
    /// # Errors
    ///
    /// Returns an error when the matrix dimensions are inconsistent, when a
    /// nonzero cross-component velocity entry is present, or when a component
    /// block cannot be factored by the provider sparse LU implementation.
    pub fn new(matrix: &SparseMatrix<T>, n_velocity: usize, n_pressure: usize) -> Result<Self> {
        Self::new_with_symbols(matrix, n_velocity, n_pressure, None)
            .map(|(preconditioner, _)| preconditioner)
    }

    /// Construct while reusing symbolic factorization for an unchanged mesh.
    pub(crate) fn new_with_cache(
        matrix: &SparseMatrix<T>,
        n_velocity: usize,
        n_pressure: usize,
        cache: &mut Option<ComponentBlockPattern>,
    ) -> Result<Self> {
        let cache_matches = cache
            .as_ref()
            .is_some_and(|pattern| pattern.matches(matrix, n_velocity, n_pressure));
        let (preconditioner, symbols) = {
            let cached_symbols = cache
                .as_ref()
                .filter(|_| cache_matches)
                .map(|pattern| pattern.symbols.as_slice());
            Self::new_with_symbols(matrix, n_velocity, n_pressure, cached_symbols)?
        };
        if !cache_matches {
            *cache = Some(ComponentBlockPattern {
                component_size: n_velocity / VELOCITY_COMPONENTS,
                n_velocity,
                n_pressure,
                row_ptr: matrix.row_ptr().to_vec(),
                col_indices: matrix.col_indices().to_vec(),
                symbols,
            });
        }
        Ok(preconditioner)
    }

    fn new_with_symbols(
        matrix: &SparseMatrix<T>,
        n_velocity: usize,
        n_pressure: usize,
        cached_symbols: Option<&[SymbolicLu]>,
    ) -> Result<(Self, Vec<SymbolicLu>)> {
        if matrix.nrows() != n_velocity + n_pressure {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Matrix size mismatch: {} != {} + {}",
                matrix.nrows(),
                n_velocity,
                n_pressure
            )));
        }
        if n_velocity == 0 || !n_velocity.is_multiple_of(VELOCITY_COMPONENTS) {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Velocity DOF count {n_velocity} is not divisible by {VELOCITY_COMPONENTS}"
            )));
        }

        let component_size = n_velocity / VELOCITY_COMPONENTS;
        let mut momentum_blocks = Vec::with_capacity(VELOCITY_COMPONENTS);
        let mut symbols = Vec::with_capacity(VELOCITY_COMPONENTS);

        for component in 0..VELOCITY_COMPONENTS {
            let offset = component * component_size;
            let mut block = SparseMatrixBuilder::new(component_size, component_size);
            for local_row in 0..component_size {
                let global_row = offset + local_row;
                let row = matrix.row(global_row);
                for (&global_col, &value) in row.col_indices().iter().zip(row.values()) {
                    if global_col >= n_velocity {
                        continue;
                    }
                    let column_component = global_col / component_size;
                    if column_component != component {
                        if NumericElement::abs(value) > diagonal_epsilon() {
                            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                                "velocity component coupling ({global_row}, {global_col}) is not supported by component-block preconditioning"
                            )));
                        }
                        continue;
                    }
                    block.add_entry(local_row, global_col - offset, value)?;
                }
            }

            let solver = SparseLuSolver {
                max_size: component_size,
                ..SparseLuSolver::default()
            };
            let block = stabilize_preconditioner_diagonal(&block.build()?, solver.pivot_tolerance)?;
            let symbol = cached_symbols
                .and_then(|symbols| symbols.get(component))
                .cloned()
                .unwrap_or_else(|| factor_symbolic(&CscMatrix::from_csr(&block)));
            let factor = solver.factor_sparse_with_symbolic(&block, &symbol)?;
            symbols.push(symbol);
            momentum_blocks.push(factor);
        }

        let simple = SimplePreconditioner::new(matrix, n_velocity, n_pressure)?;
        let mut pressure_builder = SparseMatrixBuilder::new(n_pressure, n_pressure);
        for pressure_row in 0..n_pressure {
            let row = matrix.row(n_velocity + pressure_row);
            for (&column, &value) in row.col_indices().iter().zip(row.values()) {
                if let Some(pressure_column) = column
                    .checked_sub(n_velocity)
                    .filter(|&index| index < n_pressure)
                {
                    pressure_builder.add_entry(pressure_row, pressure_column, value)?;
                }
            }
            for &(velocity, divergence_value) in &simple.divergence_rows[pressure_row] {
                let momentum_inverse = simple.momentum_inv.diag_inv[velocity];
                for pressure_column in 0..n_pressure {
                    for &(gradient_velocity, gradient_value) in
                        &simple.gradient_columns[pressure_column]
                    {
                        if gradient_velocity == velocity {
                            pressure_builder.add_entry(
                                pressure_row,
                                pressure_column,
                                -divergence_value * momentum_inverse * gradient_value,
                            )?;
                        }
                    }
                }
            }
        }
        let pressure_matrix = stabilize_preconditioner_diagonal(
            &pressure_builder.build()?,
            SparseLuSolver::default().pivot_tolerance,
        )?;
        // The pressure Schur approximation has a structural nullspace for
        // this reduced boundary topology. Do not spend a full numeric sparse
        // factorization discovering that fact on every Picard iteration; use
        // signed row-sum lumping as the explicit pressure recovery.
        let pressure_block =
            PressureBlock::LumpedDiagonal(lumped_pressure_preconditioner(&pressure_matrix));

        Ok((
            Self {
                momentum_blocks,
                pressure_block,
                simple,
                component_size,
                n_velocity,
                n_pressure,
            },
            symbols,
        ))
    }

    /// Apply the component-block momentum solves and SIMPLE pressure update.
    pub fn apply(&self, b: &Array1<T>) -> Result<Array1<T>> {
        let n_total = self.n_velocity + self.n_pressure;
        validate_vector_len("component-block preconditioner input", b, n_total)?;

        let mut u_star = Array1::zeros([self.n_velocity]);
        let mut block_rhs = Array1::zeros([self.component_size]);
        let mut block_solution = Array1::zeros([self.component_size]);
        for (component, factor) in self.momentum_blocks.iter().enumerate() {
            let offset = component * self.component_size;
            for local in 0..self.component_size {
                block_rhs[local] = b[offset + local];
            }
            factor.solve_into(&block_rhs.view(), &mut block_solution.view_mut())?;
            for local in 0..self.component_size {
                u_star[offset + local] = block_solution[local];
            }
        }

        let mut rhs_p = Array1::zeros([self.n_pressure]);
        for pressure in 0..self.n_pressure {
            rhs_p[pressure] = b[self.n_velocity + pressure];
            let mut divergence_u = <T as NumericElement>::ZERO;
            for &(velocity, value) in &self.simple.divergence_rows[pressure] {
                divergence_u += value * u_star[velocity];
            }
            rhs_p[pressure] -= divergence_u;
        }

        let pressure = match &self.pressure_block {
            PressureBlock::LumpedDiagonal(preconditioner) => preconditioner.apply(&rhs_p)?,
        };

        let mut gradient_pressure = Array1::zeros([self.n_velocity]);
        for pressure_index in 0..self.n_pressure {
            for &(velocity, value) in &self.simple.gradient_columns[pressure_index] {
                gradient_pressure[velocity] += value * pressure[pressure_index];
            }
        }

        for (component, factor) in self.momentum_blocks.iter().enumerate() {
            let offset = component * self.component_size;
            for local in 0..self.component_size {
                block_rhs[local] = gradient_pressure[offset + local];
            }
            factor.solve_into(&block_rhs.view(), &mut block_solution.view_mut())?;
            for local in 0..self.component_size {
                u_star[offset + local] -= block_solution[local];
            }
        }

        let mut result = Array1::zeros([n_total]);
        for velocity in 0..self.n_velocity {
            result[velocity] = u_star[velocity];
        }
        for pressure_index in 0..self.n_pressure {
            result[self.n_velocity + pressure_index] = pressure[pressure_index];
        }
        Ok(result)
    }
}

impl<T> Preconditioner<T> for ComponentBlockPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoRealScalar,
{
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> std::result::Result<(), LetoError> {
        validate_vector_len(
            "component-block preconditioner input",
            r,
            self.n_velocity + self.n_pressure,
        )?;
        validate_vector_len("component-block preconditioner output", z, vector_len(r))?;
        let result = self
            .apply(r)
            .map_err(|error| LetoError::InvalidInput(error.to_string()))?;
        for index in 0..vector_len(z) {
            z[index] = result[index];
        }
        Ok(())
    }
}

impl<T> athena_core::Preconditioner<athena_leto::LetoBackend<T>> for ComponentBlockPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoRealScalar,
{
    /// Apply component sparse-LU solves directly to Athena's borrowed vectors.
    fn apply(
        &self,
        _backend: &athena_leto::LetoBackend<T>,
        residual: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::View<'_>,
        mut output: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::ViewMut<'_>,
    ) -> std::result::Result<(), athena_leto::LetoBackendError> {
        let expected = self.n_velocity + self.n_pressure;
        if residual.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: residual.shape()[0],
                right: expected,
            });
        }
        if output.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: output.shape()[0],
                right: expected,
            });
        }

        let mut block_rhs = Array1::zeros([self.component_size]);
        let mut block_solution = Array1::zeros([self.component_size]);
        for (component, factor) in self.momentum_blocks.iter().enumerate() {
            let offset = component * self.component_size;
            for local in 0..self.component_size {
                block_rhs[local] = residual[offset + local];
            }
            factor.solve_into(&block_rhs.view(), &mut block_solution.view_mut())?;
            for local in 0..self.component_size {
                output[offset + local] = block_solution[local];
            }
        }

        let mut pressure_rhs = Array1::zeros([self.n_pressure]);
        for pressure_index in 0..self.n_pressure {
            pressure_rhs[pressure_index] = residual[self.n_velocity + pressure_index];
            for &(velocity, value) in &self.simple.divergence_rows[pressure_index] {
                pressure_rhs[pressure_index] -= value * output[velocity];
            }
        }
        let pressure = match &self.pressure_block {
            PressureBlock::LumpedDiagonal(preconditioner) => {
                let mut pressure = Array1::zeros([self.n_pressure]);
                for index in 0..self.n_pressure {
                    pressure[index] = pressure_rhs[index] * preconditioner.diag_inv[index];
                }
                pressure
            }
        };
        for pressure_index in 0..self.n_pressure {
            output[self.n_velocity + pressure_index] = pressure[pressure_index];
        }

        for (component, factor) in self.momentum_blocks.iter().enumerate() {
            let offset = component * self.component_size;
            block_rhs.fill(<T as NumericElement>::ZERO);
            for pressure_index in 0..self.n_pressure {
                let pressure = output[self.n_velocity + pressure_index];
                for &(velocity, value) in &self.simple.gradient_columns[pressure_index] {
                    if velocity / self.component_size == component {
                        block_rhs[velocity - offset] += value * pressure;
                    }
                }
            }
            factor.solve_into(&block_rhs.view(), &mut block_solution.view_mut())?;
            for local in 0..self.component_size {
                output[offset + local] -= block_solution[local];
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparseMatrixBuilder;
    use leto::Array1;

    fn saddle_point_matrix() -> SparseMatrix<f64> {
        // Create simple 4x4 saddle-point system: 2 velocity, 2 pressure DOFs.
        let mut builder = SparseMatrixBuilder::new(4, 4);

        // Momentum block (2x2) - viscous terms.
        builder
            .add_entry(0, 0, 4.0)
            .expect("invariant: valid saddle-point fixture entry");
        builder
            .add_entry(1, 1, 4.0)
            .expect("invariant: valid saddle-point fixture entry");

        // Gradient block (velocity -> pressure coupling).
        builder
            .add_entry(0, 2, 1.0)
            .expect("invariant: valid saddle-point fixture entry");
        builder
            .add_entry(1, 3, 1.0)
            .expect("invariant: valid saddle-point fixture entry");

        // Divergence block (pressure -> velocity coupling).
        builder
            .add_entry(2, 0, 1.0)
            .expect("invariant: valid saddle-point fixture entry");
        builder
            .add_entry(3, 1, 1.0)
            .expect("invariant: valid saddle-point fixture entry");

        // Pressure block diagonal stabilization.
        builder
            .add_entry(2, 2, 1.0)
            .expect("invariant: valid saddle-point fixture entry");
        builder
            .add_entry(3, 3, 1.0)
            .expect("invariant: valid saddle-point fixture entry");

        let mut rhs = Array1::zeros([4]);
        builder
            .build_with_rhs(&mut rhs)
            .expect("invariant: valid saddle-point fixture builds")
    }

    fn normalized_saddle_point_matrix() -> SparseMatrix<f64> {
        let mut builder = SparseMatrixBuilder::new(4, 4);
        builder
            .add_entry(0, 0, 4.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .add_entry(1, 1, 4.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .add_entry(0, 2, -1.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .add_entry(1, 3, -1.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .add_entry(2, 0, -1.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .add_entry(3, 1, -1.0)
            .expect("invariant: valid normalized fixture entry");
        builder
            .build()
            .expect("invariant: valid normalized fixture builds")
    }

    fn component_saddle_point_matrix() -> SparseMatrix<f64> {
        let mut builder = SparseMatrixBuilder::new(8, 8);
        for velocity in 0..6 {
            builder
                .add_entry(velocity, velocity, 4.0)
                .expect("invariant: valid component fixture entry");
        }
        builder
            .add_entry(0, 6, -1.0)
            .expect("invariant: valid component fixture entry");
        builder
            .add_entry(1, 7, -1.0)
            .expect("invariant: valid component fixture entry");
        builder
            .add_entry(6, 0, -1.0)
            .expect("invariant: valid component fixture entry");
        builder
            .add_entry(7, 1, -1.0)
            .expect("invariant: valid component fixture entry");
        builder
            .build()
            .expect("invariant: valid component fixture builds")
    }

    #[test]
    fn test_diagonal_preconditioner() {
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder
            .add_entry(0, 0, 2.0)
            .expect("invariant: valid diagonal fixture entry");
        builder
            .add_entry(1, 1, 4.0)
            .expect("invariant: valid diagonal fixture entry");
        builder
            .add_entry(2, 2, 8.0)
            .expect("invariant: valid diagonal fixture entry");

        let mut rhs = Array1::zeros([3]);
        let matrix = builder
            .build_with_rhs(&mut rhs)
            .expect("invariant: valid diagonal fixture builds");

        let precond = DiagonalPreconditioner::new(&matrix);
        let b = Array1::from_shape_vec([3], vec![2.0, 4.0, 8.0])
            .expect("invariant: diagonal fixture shape matches values");
        let x = precond
            .apply(&b)
            .expect("invariant: diagonal fixture application succeeds");

        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);
        assert!((x[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_block_diagonal_preconditioner() {
        let matrix = saddle_point_matrix();
        let precond = BlockDiagonalPreconditioner::new(&matrix, 2, 2)
            .expect("invariant: valid block preconditioner fixture");
        let b = Array1::from_shape_vec([4], vec![4.0, 4.0, 1.0, 1.0])
            .expect("invariant: block fixture shape matches values");
        let x = precond
            .apply(&b)
            .expect("invariant: block fixture application succeeds");

        // Momentum block: u = [4/4, 4/4] = [1, 1]
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);

        // Pressure block: p = [1/1, 1/1] = [1, 1]
        assert!((x[2] - 1.0).abs() < 1e-10);
        assert!((x[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn simple_preconditioner_uses_leto_arrays_for_coupled_correction() {
        let matrix = saddle_point_matrix();
        let precond = SimplePreconditioner::new(&matrix, 2, 2)
            .expect("invariant: valid SIMPLE preconditioner fixture");
        let b = Array1::from_shape_vec([4], vec![4.0, 4.0, 2.0, 2.0])
            .expect("invariant: SIMPLE fixture shape matches values");
        let x = precond
            .apply(&b)
            .expect("invariant: SIMPLE fixture application succeeds");

        // C - D diag(A)^-1 G = 3/4, yielding the exact solution.
        assert!((x[0] - 2.0 / 3.0).abs() < 1e-10);
        assert!((x[1] - 2.0 / 3.0).abs() < 1e-10);
        assert!((x[2] - 4.0 / 3.0).abs() < 1e-10);
        assert!((x[3] - 4.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn simple_preconditioner_preserves_normalized_continuity_sign() {
        let matrix = normalized_saddle_point_matrix();
        let precond = SimplePreconditioner::new(&matrix, 2, 2)
            .expect("invariant: valid normalized SIMPLE fixture");
        let b = Array1::from_shape_vec([4], vec![4.0, 4.0, 2.0, 2.0])
            .expect("invariant: normalized SIMPLE shape matches values");
        let x = precond
            .apply(&b)
            .expect("invariant: normalized SIMPLE application succeeds");

        assert!((x[0] + 2.0).abs() < 1e-10);
        assert!((x[1] + 2.0).abs() < 1e-10);
        assert!((x[2] + 12.0).abs() < 1e-10);
        assert!((x[3] + 12.0).abs() < 1e-10);
    }

    #[test]
    fn component_block_preconditioner_factors_provider_velocity_blocks() {
        let matrix = component_saddle_point_matrix();
        let precond = ComponentBlockPreconditioner::new(&matrix, 6, 2)
            .expect("invariant: valid component preconditioner fixture");
        let b = Array1::from_shape_vec([8], vec![4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0])
            .expect("invariant: component fixture shape matches values");
        let x = precond
            .apply(&b)
            .expect("invariant: component fixture application succeeds");

        assert!((x[0] + 2.0).abs() < 1e-10);
        assert!((x[1] + 2.0).abs() < 1e-10);
        assert!(x[2].abs() < 1e-10);
        assert!(x[3].abs() < 1e-10);
        assert!(x[4].abs() < 1e-10);
        assert!(x[5].abs() < 1e-10);
        assert!((x[6] + 12.0).abs() < 1e-10);
        assert!((x[7] + 12.0).abs() < 1e-10);
    }

    #[test]
    fn preconditioners_reject_mismatched_vector_lengths() {
        let matrix = saddle_point_matrix();
        let block = BlockDiagonalPreconditioner::new(&matrix, 2, 2)
            .expect("invariant: valid block preconditioner fixture");
        let simple = SimplePreconditioner::new(&matrix, 2, 2)
            .expect("invariant: valid SIMPLE preconditioner fixture");
        let wrong = Array1::zeros([3]);

        block
            .apply(&wrong)
            .expect_err("block preconditioner must reject a vector length mismatch");
        simple
            .apply(&wrong)
            .expect_err("SIMPLE preconditioner must reject a vector length mismatch");
    }
}

impl<T> athena_core::Preconditioner<athena_leto::LetoBackend<T>> for BlockDiagonalPreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoRealScalar,
{
    /// Block-diagonal apply straight over the borrowed views.
    ///
    /// The velocity and pressure blocks are contiguous index ranges, so this
    /// needs no scratch and no copy.
    fn apply(
        &self,
        _backend: &athena_leto::LetoBackend<T>,
        residual: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::View<'_>,
        mut output: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::ViewMut<'_>,
    ) -> std::result::Result<(), athena_leto::LetoBackendError> {
        let expected = self.n_velocity + self.n_pressure;
        if residual.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: residual.shape()[0],
                right: expected,
            });
        }
        if output.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: output.shape()[0],
                right: expected,
            });
        }
        for index in 0..self.n_velocity {
            output[index] = residual[index] * self.momentum_preconditioner.diag_inv[index];
        }
        for index in 0..self.n_pressure {
            let offset = self.n_velocity + index;
            output[offset] = residual[offset] * self.pressure_preconditioner.diag_inv[index];
        }
        Ok(())
    }
}

impl<T> athena_core::Preconditioner<athena_leto::LetoBackend<T>> for SimplePreconditioner<T>
where
    T: RealField + FloatElement + Copy + LetoRealScalar,
{
    /// Apply SIMPLE directly to Athena's borrowed vectors.
    ///
    /// The velocity part of `output` first stores the diagonal momentum
    /// prediction. Pressure is then formed from the borrowed divergence rows,
    /// and the same output buffer receives the pressure correction. This
    /// preserves the SIMPLE recurrence without allocating an intermediate
    /// vector on the Krylov hot path.
    fn apply(
        &self,
        _backend: &athena_leto::LetoBackend<T>,
        residual: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::View<'_>,
        mut output: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::ViewMut<'_>,
    ) -> std::result::Result<(), athena_leto::LetoBackendError> {
        let expected = self.n_velocity + self.n_pressure;
        if residual.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: residual.shape()[0],
                right: expected,
            });
        }
        if output.shape()[0] != expected {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: output.shape()[0],
                right: expected,
            });
        }

        for index in 0..self.n_velocity {
            output[index] = residual[index] * self.momentum_inv.diag_inv[index];
        }
        for pressure_index in 0..self.n_pressure {
            let mut pressure_rhs = residual[self.n_velocity + pressure_index];
            for &(velocity_index, divergence_value) in &self.divergence_rows[pressure_index] {
                pressure_rhs -= divergence_value * output[velocity_index];
            }
            let pressure_offset = self.n_velocity + pressure_index;
            output[pressure_offset] = pressure_rhs * self.schur_diag_inv[pressure_index];
        }
        for pressure_index in 0..self.n_pressure {
            let pressure = output[self.n_velocity + pressure_index];
            for &(velocity_index, gradient_value) in &self.gradient_columns[pressure_index] {
                output[velocity_index] -=
                    self.momentum_inv.diag_inv[velocity_index] * gradient_value * pressure;
            }
        }
        Ok(())
    }
}
