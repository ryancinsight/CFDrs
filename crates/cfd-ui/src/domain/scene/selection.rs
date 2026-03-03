//! Selection tracking for scene nodes.

use std::collections::HashSet;

/// The selection mode for picking operations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SelectionMode {
    /// Replace the current selection.
    Replace,
    /// Toggle the clicked item in the selection.
    Toggle,
    /// Add to the current selection.
    Add,
}

/// Tracks which scene nodes are currently selected.
#[derive(Clone, Debug, Default)]
pub struct SelectionSet {
    selected: HashSet<usize>,
}

impl SelectionSet {
    /// Create an empty selection.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Select a single node, replacing any existing selection.
    pub fn select(&mut self, idx: usize) {
        self.selected.clear();
        self.selected.insert(idx);
    }

    /// Toggle a node in the selection.
    pub fn toggle(&mut self, idx: usize) {
        if !self.selected.remove(&idx) {
            self.selected.insert(idx);
        }
    }

    /// Add a node to the selection without clearing.
    pub fn add(&mut self, idx: usize) {
        self.selected.insert(idx);
    }

    /// Clear the entire selection.
    pub fn clear(&mut self) {
        self.selected.clear();
    }

    /// Check if a node is selected.
    #[must_use]
    pub fn is_selected(&self, idx: usize) -> bool {
        self.selected.contains(&idx)
    }

    /// Number of selected items.
    #[must_use]
    pub fn count(&self) -> usize {
        self.selected.len()
    }

    /// Iterate over selected node indices.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.selected.iter().copied()
    }

    /// Apply a selection operation based on mode.
    pub fn apply(&mut self, idx: usize, mode: SelectionMode) {
        match mode {
            SelectionMode::Replace => self.select(idx),
            SelectionMode::Toggle => self.toggle(idx),
            SelectionMode::Add => self.add(idx),
        }
    }

    /// Get the single selected node, if exactly one is selected.
    #[must_use]
    pub fn single(&self) -> Option<usize> {
        if self.selected.len() == 1 {
            self.selected.iter().next().copied()
        } else {
            None
        }
    }
}
