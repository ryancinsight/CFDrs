//! Mesh connectivity information.

/// Mesh connectivity structure
pub struct Connectivity;

impl Default for Connectivity {
    fn default() -> Self {
        Self::new()
    }
}

impl Connectivity {
    /// Create a new connectivity structure
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}
