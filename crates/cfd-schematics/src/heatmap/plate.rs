/// Scale factor: 1 mm → 6 px.
const SCALE: f64 = 6.0;

/// Left/right margin [px].
pub(super) const MARGIN_X: f64 = 62.0;
/// Top margin [px].
pub(super) const MARGIN_Y: f64 = 70.0;

/// SBS plate width \[mm].
pub(super) const PLATE_W_MM: f64 = 127.76;
/// SBS plate height \[mm].
pub(super) const PLATE_H_MM: f64 = 85.47;

/// Well pitch \[mm].
pub(super) const PITCH: f64 = 9.0;
/// First well (A1) centre X \[mm].
pub(super) const WELL_A1_X: f64 = 14.38;
/// First well (A1) centre Y \[mm].
pub(super) const WELL_A1_Y: f64 = 11.24;
/// Well drawing radius \[mm].
pub(super) const WELL_R: f64 = 3.5;

/// Treatment-zone centre-to-centre span \[mm] across 6 wells (5 pitches).
pub(super) const ZONE_CENTER_SPAN_MM: f64 = 45.0;
/// Treatment zone first column index (0-based) — column 4 in 1-indexed (3 in 0-indexed).
pub(super) const ZONE_COL_START: usize = 3;
/// Treatment zone first row index (0-based) — row B in A–H labelling (1 in 0-indexed).
pub(super) const ZONE_ROW_START: usize = 1;
/// Number of wells in each axis of the treatment zone.
pub(super) const ZONE_WELLS: usize = 6;
/// Full highlighted treatment-zone envelope \[mm] including a half-pitch border
/// around the first and last well centres.
pub(super) const ZONE_ENVELOPE_MM: f64 = PITCH * ZONE_WELLS as f64;

#[inline]
pub(super) fn mm_to_px(mm: f64) -> f64 {
    mm * SCALE
}

#[inline]
pub(super) fn px_x(mm: f64) -> f64 {
    MARGIN_X + mm_to_px(mm)
}

#[inline]
pub(super) fn px_y(mm: f64) -> f64 {
    MARGIN_Y + mm_to_px(mm)
}

