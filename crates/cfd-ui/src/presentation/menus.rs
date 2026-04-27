//! Application menu definitions.

/// Top-level menu structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MenuCategory {
    File,
    Edit,
    Mesh,
    Simulation,
    View,
    Help,
}

impl MenuCategory {
    /// Display label for this menu category.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::File => "File",
            Self::Edit => "Edit",
            Self::Mesh => "Mesh",
            Self::Simulation => "Simulation",
            Self::View => "View",
            Self::Help => "Help",
        }
    }

    /// All menu categories in order.
    #[must_use]
    pub fn all() -> &'static [Self] {
        &[
            Self::File,
            Self::Edit,
            Self::Mesh,
            Self::Simulation,
            Self::View,
            Self::Help,
        ]
    }
}
