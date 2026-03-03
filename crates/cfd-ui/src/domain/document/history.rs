//! Command history for undo/redo operations.

use crate::domain::document::project::ProjectDocument;

/// A command that can be executed, undone, and described.
pub trait UndoableCommand: Send + Sync {
    /// Execute the command, mutating the document.
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()>;

    /// Reverse the command, restoring the previous state.
    fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()>;

    /// Human-readable description for the undo/redo menu.
    fn description(&self) -> &str;
}

/// Tracks executed commands for undo/redo support.
pub struct CommandHistory {
    undo_stack: Vec<Box<dyn UndoableCommand>>,
    redo_stack: Vec<Box<dyn UndoableCommand>>,
}

impl CommandHistory {
    /// Create an empty history.
    #[must_use]
    pub fn new() -> Self {
        Self {
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    /// Execute a command and push it onto the undo stack.
    ///
    /// Clears the redo stack since a new action invalidates the redo chain.
    pub fn execute(
        &mut self,
        mut cmd: Box<dyn UndoableCommand>,
        doc: &mut ProjectDocument,
    ) -> anyhow::Result<()> {
        cmd.execute(doc)?;
        self.undo_stack.push(cmd);
        self.redo_stack.clear();
        Ok(())
    }

    /// Undo the most recent command.
    pub fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<bool> {
        if let Some(mut cmd) = self.undo_stack.pop() {
            cmd.undo(doc)?;
            self.redo_stack.push(cmd);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Redo the most recently undone command.
    pub fn redo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<bool> {
        if let Some(mut cmd) = self.redo_stack.pop() {
            cmd.execute(doc)?;
            self.undo_stack.push(cmd);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Whether there is a command that can be undone.
    #[must_use]
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// Whether there is a command that can be redone.
    #[must_use]
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Description of the next undo command (for menu label).
    #[must_use]
    pub fn undo_description(&self) -> Option<&str> {
        self.undo_stack.last().map(|c| c.description())
    }

    /// Description of the next redo command (for menu label).
    #[must_use]
    pub fn redo_description(&self) -> Option<&str> {
        self.redo_stack.last().map(|c| c.description())
    }
}

impl Default for CommandHistory {
    fn default() -> Self {
        Self::new()
    }
}
