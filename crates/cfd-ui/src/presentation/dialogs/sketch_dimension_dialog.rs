//! Sketch dimension dialog — enter dimension value and toggle driving/driven.

/// State for the dimension input dialog.
#[derive(Clone, Debug)]
pub struct DimensionDialogState {
    /// Current input value (in model units).
    pub value: f64,
    /// Whether the dimension is driving (constrains geometry) or driven (read-only).
    pub driving: bool,
    /// Label shown in the dialog title.
    pub label: String,
    /// Whether the dialog is currently visible.
    pub visible: bool,
}

impl DimensionDialogState {
    /// Create a new dialog state for a distance dimension.
    #[must_use]
    pub fn distance(initial_value: f64) -> Self {
        Self {
            value: initial_value,
            driving: true,
            label: "Distance".to_owned(),
            visible: true,
        }
    }

    /// Create a new dialog state for an angle dimension.
    #[must_use]
    pub fn angle(initial_value_deg: f64) -> Self {
        Self {
            value: initial_value_deg,
            driving: true,
            label: "Angle".to_owned(),
            visible: true,
        }
    }

    /// Create a new dialog state for a radius dimension.
    #[must_use]
    pub fn radius(initial_value: f64) -> Self {
        Self {
            value: initial_value,
            driving: true,
            label: "Radius".to_owned(),
            visible: true,
        }
    }

    /// Dismiss the dialog.
    pub fn dismiss(&mut self) {
        self.visible = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_dialog_is_driving_by_default() {
        let d = DimensionDialogState::distance(5.0);
        assert!(d.driving);
        assert!(d.visible);
        assert_eq!(d.value, 5.0);
    }

    #[test]
    fn dismiss_hides_dialog() {
        let mut d = DimensionDialogState::radius(1.0);
        d.dismiss();
        assert!(!d.visible);
    }
}
