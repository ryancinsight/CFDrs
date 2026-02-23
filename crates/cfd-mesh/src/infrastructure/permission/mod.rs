//! GhostCell-inspired permission system.
//!
//! Separates **data** from **access rights** at compile time. Mesh data lives
//! inside `GhostCell<'brand, T>` wrappers. Reading or writing requires a
//! `GhostToken<'brand>` whose lifetime brand must match. Since the brand is
//! created via a closure (`GhostToken::new(|token| { ... })`), it is
//! impossible to forge a token from an unrelated scope.
//!
//! This gives us:
//! - Zero-cost abstraction (no runtime checks)
//! - Single-writer OR multiple-reader at compile time
//! - Safe interior mutability without `RefCell`/`Mutex` overhead
//!
//! ## Usage
//!
//! ```rust,ignore
//! use cfd_mesh::permission::{GhostToken, GhostCell};
//!
//! GhostToken::new(|mut token| {
//!     let cell = GhostCell::new(42u32);
//!     assert_eq!(*cell.borrow(&token), 42);
//!     *cell.borrow_mut(&mut token) = 99;
//!     assert_eq!(*cell.borrow(&token), 99);
//! });
//! ```

pub mod token;
pub mod cell;
pub mod arena;

pub use token::GhostToken;
pub use cell::GhostCell;
pub use arena::PermissionedArena;
