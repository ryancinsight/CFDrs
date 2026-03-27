//! Selection module — tracks which scene entities are currently selected.

pub mod granularity;
pub mod set;
pub mod target;

pub use granularity::SelectionGranularity;
pub use set::{SelectionMode, SelectionSet};
pub use target::SelectionTarget;
