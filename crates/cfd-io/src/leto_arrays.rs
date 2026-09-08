use leto::Array2;

pub(crate) fn try_for_each_row_major<T, E>(
    array: &Array2<T>,
    mut visit: impl FnMut(T) -> Result<(), E>,
) -> Result<(), E>
where
    T: Copy,
{
    let [rows, cols] = array.shape();
    for row in 0..rows {
        for col in 0..cols {
            let value = *array
                .get([row, col])
                .expect("invariant: generated row-major Leto index is in bounds");
            visit(value)?;
        }
    }
    Ok(())
}

pub(crate) fn row_major_values<T>(array: &Array2<T>) -> Vec<T>
where
    T: Copy,
{
    let mut values = Vec::with_capacity(array.size());
    try_for_each_row_major(array, |value| {
        values.push(value);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("invariant: infallible row-major collection cannot fail");
    values
}

pub(crate) fn all_row_major<T>(array: &Array2<T>, mut predicate: impl FnMut(T) -> bool) -> bool
where
    T: Copy,
{
    let [rows, cols] = array.shape();
    for row in 0..rows {
        for col in 0..cols {
            let value = *array
                .get([row, col])
                .expect("invariant: generated row-major Leto index is in bounds");
            if !predicate(value) {
                return false;
            }
        }
    }
    true
}
