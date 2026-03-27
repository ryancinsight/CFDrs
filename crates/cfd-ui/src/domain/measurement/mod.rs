//! Measurement domain module — types for interactive measurement tools.

pub mod query;
pub mod result;

pub use query::{MeasureToolState, MeasurementQuery};
pub use result::{MeasurementId, MeasurementKind, MeasurementResult, MeasurementStore};
