use crate::error::Result;
use cfd_core::error::Error;
use leto::{Array1, Array2};
use leto_ops::MatrixSolve;

#[cfg(test)]
pub(crate) fn vector_from_vec(values: Vec<f64>) -> Array1<f64> {
    Array1::from_shape_vec([values.len()], values)
        .expect("invariant: vector shape matches element count")
}

pub(crate) fn vector_from_element(len: usize, value: f64) -> Array1<f64> {
    Array1::from_elem([len], value)
}

pub(crate) fn vector_zeros(len: usize) -> Array1<f64> {
    Array1::zeros([len])
}

/// Overwrite every element of `vector` with `value`.
///
/// Reusing an accumulator across loop iterations needs its previous state
/// cleared; this is what the discarded fresh allocation used to provide.
pub(crate) fn vector_fill(vector: &mut Array1<f64>, value: f64) {
    for index in 0..vector.shape()[0] {
        vector[index] = value;
    }
}

pub(crate) fn vector_len(vector: &Array1<f64>) -> usize {
    vector.shape()[0]
}

pub(crate) fn vector_sum(vector: &Array1<f64>) -> f64 {
    vector.iter().copied().sum()
}

pub(crate) fn vector_norm(vector: &Array1<f64>) -> f64 {
    vector
        .iter()
        .map(|&value| value * value)
        .sum::<f64>()
        .sqrt()
}

pub(crate) fn vector_amax(vector: &Array1<f64>) -> f64 {
    vector.iter().map(|&value| value.abs()).fold(0.0, f64::max)
}

pub(crate) fn vector_dot(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector dot operands must have equal length"
    );
    lhs.iter().zip(rhs.iter()).map(|(&l, &r)| l * r).sum()
}

pub(crate) fn vector_add(lhs: &Array1<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector add operands must have equal length"
    );
    Array1::from_shape_fn(lhs.shape(), |idx| lhs[idx] + rhs[idx])
}

pub(crate) fn vector_sub(lhs: &Array1<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: vector sub operands must have equal length"
    );
    Array1::from_shape_fn(lhs.shape(), |idx| lhs[idx] - rhs[idx])
}

pub(crate) fn vector_scale(vector: &Array1<f64>, scale: f64) -> Array1<f64> {
    Array1::from_shape_fn(vector.shape(), |idx| vector[idx] * scale)
}

pub(crate) fn matrix_zeros(rows: usize, cols: usize) -> Array2<f64> {
    Array2::zeros([rows, cols])
}

#[cfg(test)]
pub(crate) fn matrix_from_vec(rows: usize, cols: usize, values: Vec<f64>) -> Array2<f64> {
    Array2::from_shape_vec([rows, cols], values)
        .expect("invariant: matrix shape matches element count")
}

#[cfg(test)]
pub(crate) fn matrix_from_element(rows: usize, cols: usize, value: f64) -> Array2<f64> {
    Array2::from_elem([rows, cols], value)
}

pub(crate) fn matrix_rows(matrix: &Array2<f64>) -> usize {
    matrix.shape()[0]
}

pub(crate) fn matrix_cols(matrix: &Array2<f64>) -> usize {
    matrix.shape()[1]
}

pub(crate) fn matrix_len(matrix: &Array2<f64>) -> usize {
    let [rows, cols] = matrix.shape();
    rows * cols
}

pub(crate) fn matrix_norm(matrix: &Array2<f64>) -> f64 {
    matrix
        .iter()
        .map(|&value| value * value)
        .sum::<f64>()
        .sqrt()
}

pub(crate) fn matrix_neg(matrix: &Array2<f64>) -> Array2<f64> {
    Array2::from_shape_fn(matrix.shape(), |idx| -matrix[idx])
}

pub(crate) fn matrix_add(lhs: &Array2<f64>, rhs: &Array2<f64>) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix add operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] + rhs[idx])
}

pub(crate) fn matrix_sub(lhs: &Array2<f64>, rhs: &Array2<f64>) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix sub operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] - rhs[idx])
}

pub(crate) fn matrix_scale(matrix: &Array2<f64>, scale: f64) -> Array2<f64> {
    Array2::from_shape_fn(matrix.shape(), |idx| matrix[idx] * scale)
}

pub(crate) fn matrix_add_scaled(lhs: &Array2<f64>, rhs: &Array2<f64>, scale: f64) -> Array2<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: matrix add-scaled operands must have equal shape"
    );
    Array2::from_shape_fn(lhs.shape(), |idx| lhs[idx] + scale * rhs[idx])
}

pub(crate) fn matrix_add_assign_scaled(target: &mut Array2<f64>, source: &Array2<f64>, scale: f64) {
    assert_eq!(
        target.shape(),
        source.shape(),
        "invariant: matrix add-assign operands must have equal shape"
    );
    let [rows, cols] = target.shape();
    for row in 0..rows {
        for col in 0..cols {
            target[[row, col]] += scale * source[[row, col]];
        }
    }
}

pub(crate) fn matrix_identity(size: usize) -> Array2<f64> {
    Array2::from_shape_fn(
        [size, size],
        |[row, col]| if row == col { 1.0 } else { 0.0 },
    )
}

pub(crate) fn matrix_flatten(matrix: &Array2<f64>) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    Array1::from_shape_fn([rows * cols], |[idx]| matrix[[idx / cols, idx % cols]])
}

pub(crate) fn matrix_from_flat(rows: usize, cols: usize, values: &Array1<f64>) -> Array2<f64> {
    assert_eq!(
        values.shape(),
        [rows * cols],
        "invariant: flattened vector length must match matrix shape"
    );
    Array2::from_shape_fn([rows, cols], |[row, col]| values[row * cols + col])
}

pub(crate) fn matrix_solve(matrix: &Array2<f64>, rhs: &Array1<f64>) -> Result<Array1<f64>> {
    matrix.solve(&rhs.view()).map_err(|err| {
        Error::Solver(format!(
            "Leto dense solve failed for {}x{} matrix and RHS length {}: {err}",
            matrix_rows(matrix),
            matrix_cols(matrix),
            vector_len(rhs)
        ))
    })
}

pub(crate) fn matrix_transpose_vector_mul(
    matrix: &Array2<f64>,
    vector: &Array1<f64>,
) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert_eq!(
        vector.shape(),
        [rows],
        "invariant: transposed matrix-vector dimensions must match"
    );
    Array1::from_shape_fn([cols], |[col]| {
        (0..rows).map(|row| matrix[[row, col]] * vector[row]).sum()
    })
}

pub(crate) fn row_vector(matrix: &Array2<f64>, row: usize) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert!(row < rows, "invariant: requested row is in bounds");
    Array1::from_shape_fn([cols], |[col]| matrix[[row, col]])
}

pub(crate) fn set_row(matrix: &mut Array2<f64>, row: usize, values: &Array1<f64>) {
    let [rows, cols] = matrix.shape();
    assert!(row < rows, "invariant: target row is in bounds");
    assert_eq!(
        values.shape(),
        [cols],
        "invariant: row value length must match matrix columns"
    );
    for col in 0..cols {
        matrix[[row, col]] = values[col];
    }
}

pub(crate) fn set_column(matrix: &mut Array2<f64>, col: usize, values: &Array1<f64>) {
    let [rows, cols] = matrix.shape();
    assert!(col < cols, "invariant: target column is in bounds");
    assert_eq!(
        values.shape(),
        [rows],
        "invariant: column value length must match matrix rows"
    );
    for row in 0..rows {
        matrix[[row, col]] = values[row];
    }
}

/// Accumulate `scale` times `matrix`'s column `col` into `target`.
///
/// This replaced a pair that copied the column first: reading it in place
/// performs the same additions over the same values in the same order, so the
/// result is bitwise identical, and the callers -- inside basis and
/// quadrature loops -- stop allocating a vector per call to read values they
/// drop immediately.
pub(crate) fn vector_add_assign_scaled_column(
    target: &mut Array1<f64>,
    matrix: &Array2<f64>,
    col: usize,
    scale: f64,
) {
    let [rows, cols] = matrix.shape();
    assert!(col < cols, "invariant: requested column is in bounds");
    assert_eq!(
        target.shape()[0],
        rows,
        "invariant: add-assign operands must have equal length"
    );
    for row in 0..rows {
        target[row] += scale * matrix[[row, col]];
    }
}

pub(crate) fn column_vector(matrix: &Array2<f64>, col: usize) -> Array1<f64> {
    let [rows, cols] = matrix.shape();
    assert!(col < cols, "invariant: requested column is in bounds");
    Array1::from_shape_fn([rows], |[row]| matrix[[row, col]])
}
