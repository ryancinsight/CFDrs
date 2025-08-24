//! Zero-copy iterator utilities for CFD operations
//!
//! This module provides efficient iterator combinators optimized for CFD computations,
//! strictly following zero-copy principles with proper borrowing.

mod norms;
mod statistics;
mod windows;
mod stencils;
mod parallel;

pub use norms::NormIteratorExt;
pub use statistics::StatisticsIteratorExt;
pub use windows::{WindowIterator, StridedWindowIterator};
pub use stencils::{StencilIterator, StencilPattern};
pub use parallel::ParallelIteratorExt;

/// Extension trait combining all iterator operations
pub trait MathIteratorExt: NormIteratorExt + StatisticsIteratorExt {}

impl<T: Iterator> MathIteratorExt for T {}