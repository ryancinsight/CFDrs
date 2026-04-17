//! Clipping commands — undoable operations on the clip plane set.
//!
//! These commands operate on `ClipPlaneSet` stored in the workspace,
//! not in the `ProjectDocument`. Since `UndoableCommand` takes
//! `&mut ProjectDocument`, the clip plane set is stored separately
//! and these commands are applied directly (not through the history).
//! This is consistent with `SolidWorks` where section planes are view
//! state, not model state.

use crate::domain::clipping::{ClipPlane, ClipPlaneId, ClipPlaneSet};

/// Add a new clip plane to the set.
#[derive(Clone, Debug)]
pub struct AddClipPlaneCommand {
    plane: ClipPlane,
    assigned_id: Option<ClipPlaneId>,
}

impl AddClipPlaneCommand {
    /// Create a new add command.
    #[must_use]
    pub fn new(plane: ClipPlane) -> Self {
        Self {
            plane,
            assigned_id: None,
        }
    }

    /// Execute: add the plane to the set. Returns the assigned ID, or `None`
    /// if the set is full.
    pub fn execute(&mut self, set: &mut ClipPlaneSet) -> Option<ClipPlaneId> {
        match set.add(self.plane.clone()) {
            Ok(id) => {
                self.assigned_id = Some(id);
                Some(id)
            }
            Err(_) => None,
        }
    }

    /// Undo: remove the plane that was added.
    pub fn undo(&self, set: &mut ClipPlaneSet) {
        if let Some(id) = self.assigned_id {
            set.remove(id);
        }
    }
}

/// Remove a clip plane from the set.
#[derive(Clone, Debug)]
pub struct RemoveClipPlaneCommand {
    id: ClipPlaneId,
    removed: Option<ClipPlane>,
}

impl RemoveClipPlaneCommand {
    /// Create a new remove command.
    #[must_use]
    pub fn new(id: ClipPlaneId) -> Self {
        Self { id, removed: None }
    }

    /// Execute: remove the plane.
    pub fn execute(&mut self, set: &mut ClipPlaneSet) {
        self.removed = set.remove(self.id);
    }

    /// Undo: re-add the removed plane.
    pub fn undo(&self, set: &mut ClipPlaneSet) {
        if let Some(ref plane) = self.removed {
            let _ = set.add(plane.clone());
        }
    }
}

/// Update a clip plane's properties.
#[derive(Clone, Debug)]
pub struct UpdateClipPlaneCommand {
    id: ClipPlaneId,
    old_plane: Option<ClipPlane>,
    new_plane: ClipPlane,
}

impl UpdateClipPlaneCommand {
    /// Create a new update command.
    #[must_use]
    pub fn new(id: ClipPlaneId, new_plane: ClipPlane) -> Self {
        Self {
            id,
            old_plane: None,
            new_plane,
        }
    }

    /// Execute: replace the plane with the new version.
    pub fn execute(&mut self, set: &mut ClipPlaneSet) {
        if let Some(existing) = set.get_mut(self.id) {
            self.old_plane = Some(existing.clone());
            existing.normal = self.new_plane.normal;
            existing.offset = self.new_plane.offset;
            existing.enabled = self.new_plane.enabled;
            existing.visible = self.new_plane.visible;
            existing.name.clone_from(&self.new_plane.name);
            existing.color = self.new_plane.color;
        }
    }

    /// Undo: restore the old plane properties.
    pub fn undo(&self, set: &mut ClipPlaneSet) {
        if let Some(ref old) = self.old_plane {
            if let Some(existing) = set.get_mut(self.id) {
                existing.normal = old.normal;
                existing.offset = old.offset;
                existing.enabled = old.enabled;
                existing.visible = old.visible;
                existing.name.clone_from(&old.name);
                existing.color = old.color;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::clipping::{ClipPlaneSet, ClipPreset};

    #[test]
    fn add_and_undo() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let mut cmd = AddClipPlaneCommand::new(plane);
        let id = cmd.execute(&mut set).expect("set not full");
        assert_eq!(set.count(), 1);
        cmd.undo(&mut set);
        assert_eq!(set.count(), 0);
        assert_eq!(id.0, 0);
    }

    #[test]
    fn remove_and_undo() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0);
        let id = set.add(plane).expect("should add");
        let mut cmd = RemoveClipPlaneCommand::new(id);
        cmd.execute(&mut set);
        assert_eq!(set.count(), 0);
        cmd.undo(&mut set);
        assert_eq!(set.count(), 1);
    }
}
