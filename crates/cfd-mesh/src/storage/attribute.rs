//! Per-element attribute storage.
//!
//! Attributes (scalar fields, region tags, boundary conditions) are stored
//! separately from mesh topology. This follows SoC: the mesh topology module
//! doesn't know about attribute semantics, and attribute storage doesn't
//! know about topology.

use std::collections::HashMap;

/// A named collection of per-element scalar attributes.
///
/// Generic over the index type `I` (VertexId, FaceId, etc.) so the same
/// mechanism works for vertex attributes, face attributes, and edge attributes.
pub struct AttributeStore<I: std::hash::Hash + Eq + Copy> {
    /// Named attribute channels. Each channel maps element index → value.
    channels: HashMap<String, HashMap<I, f64>>,
}

impl<I: std::hash::Hash + Eq + Copy> AttributeStore<I> {
    /// Create an empty attribute store.
    pub fn new() -> Self {
        Self {
            channels: HashMap::new(),
        }
    }

    /// Set a scalar attribute for an element.
    pub fn set(&mut self, channel: &str, id: I, value: f64) {
        self.channels
            .entry(channel.to_string())
            .or_default()
            .insert(id, value);
    }

    /// Get a scalar attribute for an element.
    pub fn get(&self, channel: &str, id: I) -> Option<f64> {
        self.channels.get(channel).and_then(|ch| ch.get(&id).copied())
    }

    /// Check if a channel exists.
    pub fn has_channel(&self, channel: &str) -> bool {
        self.channels.contains_key(channel)
    }

    /// List all channel names.
    pub fn channel_names(&self) -> Vec<&str> {
        self.channels.keys().map(|s| s.as_str()).collect()
    }

    /// Remove a channel.
    pub fn remove_channel(&mut self, channel: &str) {
        self.channels.remove(channel);
    }

    /// Clear all attributes.
    pub fn clear(&mut self) {
        self.channels.clear();
    }

    /// Number of channels.
    pub fn channel_count(&self) -> usize {
        self.channels.len()
    }
}

impl<I: std::hash::Hash + Eq + Copy> Default for AttributeStore<I> {
    fn default() -> Self {
        Self::new()
    }
}
