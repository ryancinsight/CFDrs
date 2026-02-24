//! Permissioned arena — a growable pool of `GhostCell`-wrapped elements.
//!
//! This is the building block for vertex/face/edge pools where the data
//! is stored contiguously but access is gated by a branded token.

use super::cell::GhostCell;
use super::token::GhostToken;

/// A contiguous arena of `GhostCell`-wrapped values.
///
/// Elements are appended and accessed by index. The arena itself can be
/// shared, but reading/writing individual elements requires a
/// `GhostToken<'brand>`.
pub struct PermissionedArena<'brand, T> {
    elements: Vec<GhostCell<'brand, T>>,
}

impl<'brand, T> PermissionedArena<'brand, T> {
    /// Create an empty arena.
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
        }
    }

    /// Create an arena with pre-allocated capacity.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            elements: Vec::with_capacity(cap),
        }
    }

    /// Number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// Is the arena empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Push a new element, returning its index.
    pub fn push(&mut self, value: T) -> usize {
        let idx = self.elements.len();
        self.elements.push(GhostCell::new(value));
        idx
    }

    /// Borrow element at `index` immutably.
    ///
    /// # Panics
    /// Panics if `index >= len()`.
    #[inline]
    pub fn get<'a>(&'a self, index: usize, token: &'a GhostToken<'brand>) -> &'a T {
        self.elements[index].borrow(token)
    }

    /// Borrow element at `index` mutably.
    ///
    /// # Panics
    /// Panics if `index >= len()`.
    #[inline]
    pub fn get_mut<'a>(&'a self, index: usize, token: &'a mut GhostToken<'brand>) -> &'a mut T {
        self.elements[index].borrow_mut(token)
    }

    /// Try to borrow element at `index` immutably.
    #[inline]
    pub fn try_get<'a>(&'a self, index: usize, token: &'a GhostToken<'brand>) -> Option<&'a T> {
        self.elements.get(index).map(|c| c.borrow(token))
    }

    /// Iterate over all elements immutably.
    pub fn iter<'a>(
        &'a self,
        token: &'a GhostToken<'brand>,
    ) -> impl Iterator<Item = &'a T> + use<'a, 'brand, T> {
        self.elements.iter().map(move |c| c.borrow(token))
    }

    /// Iterate with indices.
    pub fn iter_enumerated<'a>(
        &'a self,
        token: &'a GhostToken<'brand>,
    ) -> impl Iterator<Item = (usize, &'a T)> + use<'a, 'brand, T> {
        self.elements
            .iter()
            .enumerate()
            .map(move |(i, c)| (i, c.borrow(token)))
    }

    /// Clear all elements.
    pub fn clear(&mut self) {
        self.elements.clear();
    }

    /// Reserve additional capacity.
    pub fn reserve(&mut self, additional: usize) {
        self.elements.reserve(additional);
    }
}

impl<'brand, T> Default for PermissionedArena<'brand, T> {
    fn default() -> Self {
        Self::new()
    }
}
