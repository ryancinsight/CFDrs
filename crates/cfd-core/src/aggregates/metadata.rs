//! Simulation metadata

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

/// Simulation metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationMetadata {
    /// Simulation name
    pub name: String,
    /// Simulation description
    pub description: String,
    /// Creation timestamp
    pub created_at: DateTime<Utc>,
    /// Last modification timestamp
    pub modified_at: DateTime<Utc>,
    /// Version identifier
    pub version: String,
    /// Author information
    pub author: String,
    /// Custom tags for categorization
    pub tags: Vec<String>,
}

impl Default for SimulationMetadata {
    fn default() -> Self {
        let now = Utc::now();
        Self {
            name: "Unnamed Simulation".to_string(),
            description: String::new(),
            created_at: now,
            modified_at: now,
            version: "1.0.0".to_string(),
            author: String::new(),
            tags: Vec::new(),
        }
    }
}

impl SimulationMetadata {
    /// Create new metadata with name
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Set description
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = description.into();
        self.modified_at = Utc::now();
        self
    }

    /// Set author
    pub fn with_author(mut self, author: impl Into<String>) -> Self {
        self.author = author.into();
        self.modified_at = Utc::now();
        self
    }

    /// Add tag
    pub fn add_tag(&mut self, tag: impl Into<String>) {
        self.tags.push(tag.into());
        self.modified_at = Utc::now();
    }

    /// Update modification timestamp
    pub fn touch(&mut self) {
        self.modified_at = Utc::now();
    }
}
