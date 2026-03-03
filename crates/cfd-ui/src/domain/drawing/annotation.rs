//! Annotations for engineering drawings — notes, datums, surface finish symbols.

use serde::{Deserialize, Serialize};

/// A text note placed on the drawing.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NoteAnnotation {
    /// Position on the sheet in millimeters from bottom-left.
    pub position_mm: [f64; 2],
    /// The note text (may be multi-line).
    pub text: String,
    /// Font size in millimeters.
    pub font_size_mm: f64,
    /// Optional leader line endpoint (points to geometry).
    pub leader_target_mm: Option<[f64; 2]>,
}

/// A geometric dimensioning and tolerancing (GD&T) datum label.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DatumAnnotation {
    /// Position on the sheet where the datum symbol is placed.
    pub position_mm: [f64; 2],
    /// Datum letter (e.g. "A", "B", "C").
    pub datum_letter: String,
}

/// Surface finish symbol placed on the drawing.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SurfaceFinishAnnotation {
    /// Position on the sheet.
    pub position_mm: [f64; 2],
    /// Surface roughness Ra value in micrometers.
    pub roughness_ra_um: f64,
    /// Machining requirement.
    pub machining: MachiningRequirement,
}

/// Whether machining is required for a surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum MachiningRequirement {
    /// Machining required.
    Required,
    /// No machining required (as-cast, as-forged, etc.).
    NotRequired,
    /// Any process acceptable.
    Any,
}

/// Union of all annotation types.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Annotation {
    /// A text note with optional leader.
    Note(NoteAnnotation),
    /// A datum label.
    Datum(DatumAnnotation),
    /// A surface finish symbol.
    SurfaceFinish(SurfaceFinishAnnotation),
}
