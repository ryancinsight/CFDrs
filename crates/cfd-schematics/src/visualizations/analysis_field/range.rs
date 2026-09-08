//! Validated scalar range reduction.

use crate::error::{VisualizationError, VisualizationResult};

/// Precomputed finite range for constant-time normalization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct ScalarRange {
    min: f64,
    max: f64,
}

impl Default for ScalarRange {
    fn default() -> Self {
        Self { min: 0.0, max: 1.0 }
    }
}

impl ScalarRange {
    pub(super) fn from_values<'a>(
        parameter: &'static str,
        values: impl IntoIterator<Item = &'a f64>,
    ) -> VisualizationResult<Self> {
        let mut range: Option<Self> = None;

        for value in values {
            if !value.is_finite() {
                return Err(VisualizationError::InvalidParameters {
                    parameter: parameter.to_string(),
                    value: value.to_string(),
                    constraint: "all scalar field values must be finite".to_string(),
                });
            }

            range = Some(match range {
                None => Self {
                    min: *value,
                    max: *value,
                },
                Some(current) => Self {
                    min: current.min.min(*value),
                    max: current.max.max(*value),
                },
            });
        }

        Ok(range.unwrap_or_default())
    }

    pub(super) const fn endpoints(self) -> (f64, f64) {
        (self.min, self.max)
    }

    #[expect(
        clippy::float_cmp,
        reason = "exact equality detects a zero-width range; a tolerance would collapse real variation"
    )]
    pub(super) fn normalize(self, value: f64) -> f64 {
        if self.min == self.max {
            0.5
        } else {
            (value - self.min) / (self.max - self.min)
        }
    }
}
