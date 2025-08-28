//! Zero-copy iterator utilities for CFD operations
//!
//! This module provides efficient iterator combinators for CFD computations,
//! strictly following zero-copy principles with proper borrowing.

mod norms;
mod parallel;
mod statistics;
mod stencils;
mod windows;

pub use norms::NormIteratorExt;
pub use parallel::ParallelIteratorExt;
pub use statistics::StatisticsIteratorExt;
pub use stencils::{StencilIterator, StencilPattern};
pub use windows::{StridedWindowIterator, WindowIterator};

/// Extension trait combining all iterator operations
pub trait MathIteratorExt: NormIteratorExt + StatisticsIteratorExt {}

impl<T: Iterator> MathIteratorExt for T {}
