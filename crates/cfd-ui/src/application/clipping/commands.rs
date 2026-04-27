//! Clipping commands — undoable operations on the clip plane set.
//!
//! # Theorem — Slot-stable clip-plane undo
//!
//! A clipping command must restore planes at the recorded `ClipPlaneId`
//! instead of silently relocating them. If the live slot no longer matches the
//! stored snapshot, undo returns an explicit error and leaves the set unchanged.
//!
//! **Proof sketch**: `ClipPlaneSet::insert_at()` and `replace_at()` write to an
//! explicit slot only. Each command stores the plane snapshot that it added or
//! replaced and compares the live slot against that snapshot before mutating.

use crate::domain::clipping::{ClipPlane, ClipPlaneId, ClipPlaneSet, ClipPlaneSlotError};

/// Error returned by clipping commands.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ClipPlaneCommandError {
    /// The set ran out of free slots for a new plane.
    #[error("clip plane set is full")]
    SetFull,
    /// The command has no stored execution state.
    #[error("command has not been executed")]
    NotExecuted,
    /// The expected plane is absent from the set.
    #[error("clip plane slot {0:?} is missing")]
    MissingPlane(ClipPlaneId),
    /// The live slot no longer matches the stored snapshot.
    #[error("clip plane slot {0:?} changed after command execution")]
    StaleState(ClipPlaneId),
    /// The direct slot operation failed.
    #[error(transparent)]
    Slot(#[from] ClipPlaneSlotError),
}

/// Add a new clip plane to the set.
#[derive(Clone, Debug)]
pub struct AddClipPlaneCommand {
    plane: ClipPlane,
    assigned_id: Option<ClipPlaneId>,
    added_plane: Option<ClipPlane>,
}

impl AddClipPlaneCommand {
    /// Create a new add command.
    #[must_use]
    pub fn new(plane: ClipPlane) -> Self {
        Self {
            plane,
            assigned_id: None,
            added_plane: None,
        }
    }

    /// Execute: add the plane to the set. Returns the assigned ID, or an
    /// explicit error if the set is full or the recorded slot is stale.
    pub fn execute(
        &mut self,
        set: &mut ClipPlaneSet,
    ) -> Result<ClipPlaneId, ClipPlaneCommandError> {
        if let Some(id) = self.assigned_id {
            let expected = self
                .added_plane
                .as_ref()
                .ok_or(ClipPlaneCommandError::NotExecuted)?;
            match set.get(id) {
                Some(current) if current == expected => return Ok(id),
                Some(_) => return Err(ClipPlaneCommandError::StaleState(id)),
                None => {
                    set.insert_at(id, self.plane.clone())?;
                    self.added_plane = Some(
                        set.get(id)
                            .cloned()
                            .ok_or(ClipPlaneCommandError::MissingPlane(id))?,
                    );
                    return Ok(id);
                }
            }
        }

        let id = set
            .add(self.plane.clone())
            .map_err(|_| ClipPlaneCommandError::SetFull)?;
        self.assigned_id = Some(id);
        self.added_plane = Some(
            set.get(id)
                .cloned()
                .ok_or(ClipPlaneCommandError::MissingPlane(id))?,
        );
        Ok(id)
    }

    /// Undo: remove the plane that was added.
    pub fn undo(&self, set: &mut ClipPlaneSet) -> Result<(), ClipPlaneCommandError> {
        let id = self.assigned_id.ok_or(ClipPlaneCommandError::NotExecuted)?;
        let expected = self
            .added_plane
            .as_ref()
            .ok_or(ClipPlaneCommandError::NotExecuted)?;

        match set.get(id) {
            Some(current) if current == expected => {
                let removed = set
                    .remove(id)
                    .ok_or(ClipPlaneCommandError::MissingPlane(id))?;
                debug_assert_eq!(removed, *expected);
                Ok(())
            }
            None => Ok(()),
            Some(_) => Err(ClipPlaneCommandError::StaleState(id)),
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
    pub fn execute(&mut self, set: &mut ClipPlaneSet) -> Result<(), ClipPlaneCommandError> {
        if let Some(removed) = self.removed.as_ref() {
            match set.get(self.id) {
                Some(current) if current == removed => {
                    let current_removed = set
                        .remove(self.id)
                        .ok_or(ClipPlaneCommandError::MissingPlane(self.id))?;
                    debug_assert_eq!(current_removed, *removed);
                    Ok(())
                }
                None => Ok(()),
                Some(_) => Err(ClipPlaneCommandError::StaleState(self.id)),
            }
        } else {
            self.removed = Some(
                set.remove(self.id)
                    .ok_or(ClipPlaneCommandError::MissingPlane(self.id))?,
            );
            Ok(())
        }
    }

    /// Undo: re-add the removed plane.
    pub fn undo(&self, set: &mut ClipPlaneSet) -> Result<(), ClipPlaneCommandError> {
        let removed = self
            .removed
            .as_ref()
            .ok_or(ClipPlaneCommandError::NotExecuted)?;

        match set.get(self.id) {
            Some(current) if current == removed => Ok(()),
            None => {
                set.insert_at(self.id, removed.clone())?;
                Ok(())
            }
            Some(_) => Err(ClipPlaneCommandError::StaleState(self.id)),
        }
    }
}

/// Update a clip plane's properties.
#[derive(Clone, Debug)]
pub struct UpdateClipPlaneCommand {
    id: ClipPlaneId,
    old_plane: Option<ClipPlane>,
    applied_plane: Option<ClipPlane>,
    new_plane: ClipPlane,
}

impl UpdateClipPlaneCommand {
    /// Create a new update command.
    #[must_use]
    pub fn new(id: ClipPlaneId, new_plane: ClipPlane) -> Self {
        Self {
            id,
            old_plane: None,
            applied_plane: None,
            new_plane,
        }
    }

    /// Execute: replace the plane with the new version.
    pub fn execute(&mut self, set: &mut ClipPlaneSet) -> Result<(), ClipPlaneCommandError> {
        match (
            self.old_plane.as_ref(),
            self.applied_plane.as_ref(),
            set.get(self.id),
        ) {
            (Some(_old), Some(applied), Some(current)) if current == applied => Ok(()),
            (Some(old), Some(_applied), Some(current)) if current == old => {
                let previous = set.replace_at(self.id, self.new_plane.clone())?;
                debug_assert_eq!(previous, *old);
                self.applied_plane = Some(
                    set.get(self.id)
                        .cloned()
                        .ok_or(ClipPlaneCommandError::MissingPlane(self.id))?,
                );
                Ok(())
            }
            (None, None, Some(_)) => {
                let old = set.replace_at(self.id, self.new_plane.clone())?;
                self.old_plane = Some(old);
                self.applied_plane = Some(
                    set.get(self.id)
                        .cloned()
                        .ok_or(ClipPlaneCommandError::MissingPlane(self.id))?,
                );
                Ok(())
            }
            (_, _, None) => Err(ClipPlaneCommandError::MissingPlane(self.id)),
            _ => Err(ClipPlaneCommandError::StaleState(self.id)),
        }
    }

    /// Undo: restore the old plane properties.
    pub fn undo(&self, set: &mut ClipPlaneSet) -> Result<(), ClipPlaneCommandError> {
        let old = self
            .old_plane
            .as_ref()
            .ok_or(ClipPlaneCommandError::NotExecuted)?;
        let applied = self
            .applied_plane
            .as_ref()
            .ok_or(ClipPlaneCommandError::NotExecuted)?;

        match set.get(self.id) {
            Some(current) if current == old => Ok(()),
            Some(current) if current == applied => {
                let previous = set.replace_at(self.id, old.clone())?;
                debug_assert_eq!(previous, *applied);
                Ok(())
            }
            None => Err(ClipPlaneCommandError::MissingPlane(self.id)),
            Some(_) => Err(ClipPlaneCommandError::StaleState(self.id)),
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
        cmd.undo(&mut set).expect("undo should remove the plane");
        assert_eq!(set.count(), 0);
        assert_eq!(id.0, 0);
    }

    #[test]
    fn remove_and_undo() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0);
        let id = set.add(plane).expect("should add");
        let mut cmd = RemoveClipPlaneCommand::new(id);
        cmd.execute(&mut set).expect("remove should succeed");
        assert_eq!(set.count(), 0);
        cmd.undo(&mut set).expect("undo should restore the plane");
        assert_eq!(set.count(), 1);
        assert_eq!(set.get(id).expect("plane restored").id, id);
    }

    #[test]
    fn remove_undo_rejects_slot_conflict() {
        let mut set = ClipPlaneSet::new();
        let plane_a = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let plane_b = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0);
        let plane_c = ClipPlaneSet::from_preset(ClipPreset::Yz, 2.0);
        let _id_a = set.add(plane_a).expect("should add");
        let id_b = set.add(plane_b).expect("should add");
        let mut cmd = RemoveClipPlaneCommand::new(id_b);
        cmd.execute(&mut set).expect("remove should succeed");

        let occupied_id = set.add(plane_c).expect("should refill the slot");
        assert_eq!(occupied_id, id_b);

        let err = cmd
            .undo(&mut set)
            .expect_err("occupied slot should block undo");
        assert_eq!(err, ClipPlaneCommandError::StaleState(id_b));
        assert_eq!(set.get(id_b).expect("plane should remain").offset, 2.0);
    }

    #[test]
    fn update_and_undo() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let id = set.add(plane).expect("should add");
        let mut updated = ClipPlaneSet::from_preset(ClipPreset::Xy, 1.25);
        updated.enabled = false;

        let mut cmd = UpdateClipPlaneCommand::new(id, updated);
        cmd.execute(&mut set).expect("update should succeed");

        let stored = set.get(id).expect("plane stored");
        assert_eq!(stored.offset, 1.25);
        assert!(!stored.enabled);

        cmd.undo(&mut set)
            .expect("undo should restore the original plane");
        let restored = set.get(id).expect("plane restored");
        assert_eq!(restored.offset, 0.0);
        assert!(restored.enabled);
    }
}
