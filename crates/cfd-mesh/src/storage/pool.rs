//! Generic indexed pool — a contiguous `Vec<T>` with typed index access.

use std::marker::PhantomData;

/// A contiguous pool of `T` values accessed by strongly-typed index `I`.
///
/// `I` must implement `From<usize>` and provide `as_usize() -> usize`.
pub struct Pool<I, T> {
    data: Vec<T>,
    _index: PhantomData<I>,
}

impl<I, T> Pool<I, T>
where
    I: From<usize> + Copy,
{
    /// Create an empty pool.
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            _index: PhantomData,
        }
    }

    /// Create a pool with pre-allocated capacity.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            data: Vec::with_capacity(cap),
            _index: PhantomData,
        }
    }

    /// Number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Is the pool empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Insert a value, returning its typed index.
    pub fn push(&mut self, value: T) -> I {
        let idx = self.data.len();
        self.data.push(value);
        I::from(idx)
    }

    /// Access the underlying slice.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Access the underlying mutable slice.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Iterate over all elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Iterate over all elements mutably.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    /// Iterate with typed indices.
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (I, &T)> {
        self.data.iter().enumerate().map(|(i, v)| (I::from(i), v))
    }

    /// Clear all elements.
    pub fn clear(&mut self) {
        self.data.clear();
    }

    /// Reserve additional capacity.
    pub fn reserve(&mut self, additional: usize) {
        self.data.reserve(additional);
    }
}

impl<I, T> Default for Pool<I, T>
where
    I: From<usize> + Copy,
{
    fn default() -> Self {
        Self::new()
    }
}

/// Trait for typed indices that can index into a pool.
pub trait TypedIndex: From<usize> + Copy {
    /// Convert to `usize`.
    fn as_usize(self) -> usize;
}

// Implement for all our index types via a macro in the index module.
// For now, use a blanket approach — each index type's `as_usize` is available.

impl<I: TypedIndex, T> std::ops::Index<I> for Pool<I, T> {
    type Output = T;

    #[inline]
    fn index(&self, index: I) -> &T {
        &self.data[index.as_usize()]
    }
}

impl<I: TypedIndex, T> std::ops::IndexMut<I> for Pool<I, T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut T {
        &mut self.data[index.as_usize()]
    }
}
