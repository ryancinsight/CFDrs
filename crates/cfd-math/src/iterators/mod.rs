//! Zero-copy iterator utilities for CFD operations
//!
//! This module provides efficient iterator combinators for CFD computations,
//! strictly following zero-copy principles with proper borrowing.

mod statistics;
mod stencils;
mod windows;

pub use statistics::StatisticsIteratorExt;
pub use stencils::{StencilIterator, StencilPattern};
pub use windows::{StridedWindowIterator, WindowIterator};
