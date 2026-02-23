//! GhostCell — zero-cost interior mutability gated by a branded token.
//!
//! The `GhostCell<'brand, T>` wrapper stores data that can only be accessed
//! through a matching `GhostToken<'brand>`.

use std::cell::UnsafeCell;
use std::marker::PhantomData;

use super::token::GhostToken;

/// A cell whose contents are accessible only via a matching [`GhostToken`].
///
/// This is safe because:
/// - `&GhostToken` → `&T` (shared borrow of token implies shared borrow of data)
/// - `&mut GhostToken` → `&mut T` (exclusive borrow of token implies exclusive borrow of data)
///
/// The Rust borrow checker enforces that you cannot hold `&T` and `&mut T` simultaneously
/// because that would require `&GhostToken` and `&mut GhostToken` at the same time.
pub struct GhostCell<'brand, T: ?Sized> {
    _brand: PhantomData<fn(&'brand ()) -> &'brand ()>,
    value: UnsafeCell<T>,
}

// SAFETY: GhostCell is Send/Sync when T is Send/Sync because access is
// gated by the single-threaded token discipline.
// The token itself is !Send + !Sync (it uses PhantomData<fn(&'brand ()) -> &'brand ()>),
// so cross-thread access is prevented.
//
// However: for our mesh use-case, we want the *data* to be shareable across threads
// while the token stays on one thread. This is safe because the token
// prevents simultaneous mutable access.
// We mark Send+Sync only when T: Send+Sync.

// SAFETY: Access to the inner value is gated by GhostToken which ensures
// exclusive mutable access through the borrow checker. When T: Send,
// ownership transfer across threads is safe because, at any point,
// only one thread can hold &mut GhostToken.
unsafe impl<'brand, T: ?Sized + Send> Send for GhostCell<'brand, T> {}

// SAFETY: Shared references (&GhostCell) can exist on multiple threads,
// but reading requires &GhostToken (which is !Sync by construction of the
// invariant lifetime brand in PhantomData). Writing requires
// &mut GhostToken. The borrow checker prevents data races.
unsafe impl<'brand, T: ?Sized + Send + Sync> Sync for GhostCell<'brand, T> {}

impl<'brand, T> GhostCell<'brand, T> {
    /// Wrap a value in a `GhostCell`.
    #[inline]
    pub const fn new(value: T) -> Self {
        GhostCell {
            _brand: PhantomData,
            value: UnsafeCell::new(value),
        }
    }

    /// Borrow the contents immutably. Requires `&GhostToken<'brand>`.
    #[inline]
    pub fn borrow<'a>(&'a self, _token: &'a GhostToken<'brand>) -> &'a T {
        // SAFETY: The shared borrow of `_token` ensures no `&mut T` exists,
        // because obtaining `&mut T` requires `&mut GhostToken`.
        unsafe { &*self.value.get() }
    }

    /// Borrow the contents mutably. Requires `&mut GhostToken<'brand>`.
    #[inline]
    pub fn borrow_mut<'a>(&'a self, _token: &'a mut GhostToken<'brand>) -> &'a mut T {
        // SAFETY: The exclusive borrow of `_token` ensures no other `&T` or
        // `&mut T` exists, because those would require `&GhostToken` or
        // `&mut GhostToken` — but we hold `&mut GhostToken` exclusively.
        unsafe { &mut *self.value.get() }
    }

    /// Consume the cell, returning the inner value.
    #[inline]
    pub fn into_inner(self) -> T {
        self.value.into_inner()
    }
}

impl<'brand, T: Clone> GhostCell<'brand, T> {
    /// Clone the inner value. Requires `&GhostToken`.
    #[inline]
    pub fn clone_inner(&self, token: &GhostToken<'brand>) -> T {
        self.borrow(token).clone()
    }
}

impl<'brand, T: Default> Default for GhostCell<'brand, T> {
    fn default() -> Self {
        Self::new(T::default())
    }
}

impl<'brand, T: std::fmt::Debug> GhostCell<'brand, T> {
    /// Debug-format the inner value. Requires `&GhostToken`.
    pub fn debug_with(&self, token: &GhostToken<'brand>, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.borrow(token).fmt(f)
    }
}
