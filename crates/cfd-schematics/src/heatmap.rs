//! 96-well plate treatment-zone visualization.
//!
//! Generates an SVG showing a standard SBS 96-well plate (127.76 × 85.47 mm)
//! with the 6×6 SDT treatment zone highlighted and top candidate metrics
//! overlaid as colour bars inside the treatment zone.
//!
//! # SBS plate geometry (ANSI/SLAS standard)
//! - Overall: 127.76 × 85.47 mm
//! - Well pitch: 9.00 mm
//! - First well (A1) centre: x = 14.38 mm, y = 11.24 mm
//! - Grid: 12 columns × 8 rows = 96 wells
//!
//! # Treatment zone
//! - 45 × 45 mm centre-to-centre span centred at (63.88, 42.74) mm
//! - Spans wells B4–G9 (rows B–G, columns 4–9 in 1-indexed notation)
//! - The dashed overlay expands half a pitch beyond the outer wells so the
//!   full 6 × 6 well envelope is visible
//!
//! # Output
//! Pure SVG XML written to a file — no external graphical dependencies.

mod data;
mod palette;
mod plate;
mod svg;

pub use data::{CandidateZoneData, write_well_plate_diagram_svg};

#[cfg(test)]
mod tests;
