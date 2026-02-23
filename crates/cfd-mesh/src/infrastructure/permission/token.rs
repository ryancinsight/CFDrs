//! Brand-based permission token.
//!
//! The `GhostToken<'brand>` type carries a lifetime brand that ties it to a
//! specific scope. Only cells created within the same branded scope can be
//! accessed through the token.

use std::marker::PhantomData;

/// A permission token branded with lifetime `'brand`.
///
/// - `&GhostToken<'brand>` grants **shared** read access to all
///   `GhostCell<'brand, T>` values.
/// - `&mut GhostToken<'brand>` grants **exclusive** write access.
///
/// The brand ensures tokens from unrelated scopes cannot be mixed.
pub struct GhostToken<'brand> {
    _brand: PhantomData<fn(&'brand ()) -> &'brand ()>,
}

impl<'brand> GhostToken<'brand> {
    /// Create a new branded token and pass it to the closure.
    ///
    /// The brand lifetime is **invariant** — it cannot be widened or narrowed —
    /// so the token cannot escape the closure or be confused with other tokens.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// GhostToken::new(|token| {
    ///     // `token` is only valid inside this closure
    /// });
    /// ```
    pub fn new<R>(f: impl for<'new_brand> FnOnce(GhostToken<'new_brand>) -> R) -> R {
        let token = GhostToken {
            _brand: PhantomData,
        };
        f(token)
    }
}
