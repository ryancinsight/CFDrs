//! Clip plane types — Hessian normal form section planes for mesh clipping.
//!
//! # Theorem — Clip Plane Half-Space
//!
//! A clip plane in Hessian normal form `n · x + d = 0` divides R³ into two
//! half-spaces. Fragments with `n · x + d < 0` are discarded (behind the plane).
//! The signed distance from a point `p` to the plane is `n · p + d`, which is
//! positive on the visible side and negative on the clipped side.  ∎

use nalgebra::Vector3;

/// Unique identifier for a clip plane within the set (0..5).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ClipPlaneId(pub u8);

/// Error returned when a clip plane slot cannot be updated directly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ClipPlaneSlotError {
    /// The slot identifier is outside the fixed six-slot range.
    #[error("clip plane slot {0:?} is outside the 0..6 range")]
    OutOfRange(ClipPlaneId),
    /// The requested slot already contains a plane.
    #[error("clip plane slot {0:?} is already occupied")]
    Occupied(ClipPlaneId),
    /// The requested slot does not contain a plane.
    #[error("clip plane slot {0:?} is empty")]
    Empty(ClipPlaneId),
}

/// A single clip plane in Hessian normal form: `n · x + d = 0`.
///
/// Fragments on the negative side (`n · x + d < 0`) are discarded by the
/// GPU fragment shader.
#[derive(Clone, Debug, PartialEq)]
pub struct ClipPlane {
    /// Unique identifier within the clip plane set.
    pub id: ClipPlaneId,
    /// Human-readable name for the UI.
    pub name: String,
    /// Unit normal vector defining the plane orientation.
    pub normal: Vector3<f64>,
    /// Signed offset from the origin along the normal.
    pub offset: f64,
    /// Whether clipping is active for this plane.
    pub enabled: bool,
    /// Whether the plane gizmo is visible in the viewport.
    pub visible: bool,
    /// Display color for the plane gizmo (RGBA).
    pub color: [f32; 4],
}

/// Preset plane orientations aligned with coordinate axes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClipPreset {
    /// XY plane (normal along +Z).
    Xy,
    /// XZ plane (normal along +Y).
    Xz,
    /// YZ plane (normal along +X).
    Yz,
}

impl ClipPreset {
    /// Unit normal for this preset.
    #[must_use]
    fn normal(self) -> Vector3<f64> {
        match self {
            Self::Xy => Vector3::z(),
            Self::Xz => Vector3::y(),
            Self::Yz => Vector3::x(),
        }
    }

    /// Display name for this preset.
    #[must_use]
    fn label(self) -> &'static str {
        match self {
            Self::Xy => "XY Plane",
            Self::Xz => "XZ Plane",
            Self::Yz => "YZ Plane",
        }
    }

    /// Default gizmo color for this preset.
    #[must_use]
    fn color(self) -> [f32; 4] {
        match self {
            Self::Xy => [0.3, 0.3, 0.8, 0.4],
            Self::Xz => [0.3, 0.8, 0.3, 0.4],
            Self::Yz => [0.8, 0.3, 0.3, 0.4],
        }
    }
}

/// Container for up to 6 clip planes (`SolidWorks` standard).
#[derive(Clone, Debug, Default)]
pub struct ClipPlaneSet {
    planes: [Option<ClipPlane>; 6],
}

/// Error returned when the clip plane set is full (6 planes maximum).
#[derive(Clone, Debug)]
pub struct ClipPlaneSetFull;

impl ClipPlaneSet {
    /// Create an empty clip plane set.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn slot_index(id: ClipPlaneId) -> Result<usize, ClipPlaneSlotError> {
        let idx = id.0 as usize;
        if idx < 6 {
            Ok(idx)
        } else {
            Err(ClipPlaneSlotError::OutOfRange(id))
        }
    }

    /// Add a clip plane, assigning the first available slot.
    pub fn add(&mut self, mut plane: ClipPlane) -> Result<ClipPlaneId, ClipPlaneSetFull> {
        for (i, slot) in self.planes.iter_mut().enumerate() {
            if slot.is_none() {
                let id = ClipPlaneId(i as u8);
                plane.id = id;
                *slot = Some(plane);
                return Ok(id);
            }
        }
        Err(ClipPlaneSetFull)
    }

    /// Insert a clip plane at a specific slot.
    pub fn insert_at(
        &mut self,
        id: ClipPlaneId,
        mut plane: ClipPlane,
    ) -> Result<(), ClipPlaneSlotError> {
        let idx = Self::slot_index(id)?;
        let slot = &mut self.planes[idx];
        if slot.is_some() {
            return Err(ClipPlaneSlotError::Occupied(id));
        }
        plane.id = id;
        *slot = Some(plane);
        Ok(())
    }

    /// Replace the clip plane at a specific slot and return the previous plane.
    #[must_use]
    pub fn replace_at(
        &mut self,
        id: ClipPlaneId,
        mut plane: ClipPlane,
    ) -> Result<ClipPlane, ClipPlaneSlotError> {
        let idx = Self::slot_index(id)?;
        let slot = &mut self.planes[idx];
        let existing = slot.as_mut().ok_or(ClipPlaneSlotError::Empty(id))?;
        plane.id = id;
        Ok(std::mem::replace(existing, plane))
    }

    /// Remove a clip plane by ID.
    pub fn remove(&mut self, id: ClipPlaneId) -> Option<ClipPlane> {
        let idx = id.0 as usize;
        if idx < 6 {
            self.planes[idx].take()
        } else {
            None
        }
    }

    /// Access a clip plane by ID.
    #[must_use]
    pub fn get(&self, id: ClipPlaneId) -> Option<&ClipPlane> {
        self.planes.get(id.0 as usize).and_then(|s| s.as_ref())
    }

    /// Mutable access to a clip plane by ID.
    pub fn get_mut(&mut self, id: ClipPlaneId) -> Option<&mut ClipPlane> {
        self.planes.get_mut(id.0 as usize).and_then(|s| s.as_mut())
    }

    /// Iterate over all active (enabled) clip planes.
    pub fn active_planes(&self) -> impl Iterator<Item = &ClipPlane> {
        self.planes
            .iter()
            .filter_map(|s| s.as_ref().filter(|p| p.enabled))
    }

    /// Iterate over all present clip planes (enabled or not).
    pub fn all_planes(&self) -> impl Iterator<Item = &ClipPlane> {
        self.planes.iter().filter_map(|s| s.as_ref())
    }

    /// Number of active (enabled) clip planes.
    #[must_use]
    pub fn active_count(&self) -> usize {
        self.active_planes().count()
    }

    /// Number of present clip planes.
    #[must_use]
    pub fn count(&self) -> usize {
        self.planes.iter().filter(|s| s.is_some()).count()
    }

    /// Create a clip plane from a preset orientation.
    #[must_use]
    pub fn from_preset(preset: ClipPreset, offset: f64) -> ClipPlane {
        ClipPlane {
            id: ClipPlaneId(0), // assigned by add()
            name: preset.label().to_owned(),
            normal: preset.normal(),
            offset,
            enabled: true,
            visible: true,
            color: preset.color(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_and_retrieve_plane() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let id = set.add(plane).expect("should add");
        assert_eq!(id.0, 0);
        assert!(set.get(id).is_some());
        assert_eq!(set.count(), 1);
    }

    #[test]
    fn remove_plane() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0);
        let id = set.add(plane).expect("should add");
        let removed = set.remove(id);
        assert!(removed.is_some());
        assert_eq!(set.count(), 0);
    }

    #[test]
    fn set_full_returns_error() {
        let mut set = ClipPlaneSet::new();
        for _ in 0..6 {
            let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
            set.add(plane).expect("should add");
        }
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        assert!(set.add(plane).is_err());
    }

    #[test]
    fn active_count_respects_enabled() {
        let mut set = ClipPlaneSet::new();
        let mut plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        plane.enabled = false;
        set.add(plane).expect("should add");
        assert_eq!(set.active_count(), 0);
        assert_eq!(set.count(), 1);
    }

    #[test]
    fn preset_normals_are_unit() {
        use approx::assert_relative_eq;
        for preset in [ClipPreset::Xy, ClipPreset::Xz, ClipPreset::Yz] {
            assert_relative_eq!(preset.normal().norm(), 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn insert_at_preserves_requested_slot() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Yz, 0.25);
        let id = ClipPlaneId(4);
        set.insert_at(id, plane).expect("should insert");
        let stored = set.get(id).expect("plane stored");
        assert_eq!(stored.id, id);
        assert_eq!(stored.offset, 0.25);
        assert_eq!(set.count(), 1);
    }

    #[test]
    fn replace_at_returns_previous_plane() {
        let mut set = ClipPlaneSet::new();
        let first = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let second = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.25);
        let id = set.add(first.clone()).expect("should add");
        let previous = set.replace_at(id, second.clone()).expect("should replace");
        assert_eq!(previous, first);
        let stored = set.get(id).expect("plane stored");
        let mut expected = second;
        expected.id = id;
        assert_eq!(stored, &expected);
        assert_eq!(stored.id, id);
        assert_eq!(stored.offset, 1.25);
    }

    #[test]
    fn insert_at_rejects_occupied_slot() {
        let mut set = ClipPlaneSet::new();
        let id = set
            .add(ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0))
            .expect("should add");
        let err = set
            .insert_at(id, ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0))
            .expect_err("slot is occupied");
        assert_eq!(err, ClipPlaneSlotError::Occupied(id));
    }
}
