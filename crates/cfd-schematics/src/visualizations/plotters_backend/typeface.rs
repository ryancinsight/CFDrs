//! Embedded typeface for plotters text layout and rasterization.
//!
//! Plotters' `ab_glyph` font backend is pure Rust and resolves a family name
//! only against faces registered in-process, so text rendering does not depend
//! on the host's font configuration (fontconfig on Linux, DirectWrite on
//! Windows, Core Text on macOS). The crate embeds DejaVu Sans, reduced to
//! Latin, Greek, punctuation, super- and subscripts, arrows and mathematical
//! operators (`fonts/README.md` records the source and the reduction;
//! `fonts/LICENSE_DEJAVU` its license), and registers it under
//! [`EMBEDDED_FONT_FAMILY`], the family every default [`TextStyle`] names.
//!
//! A [`TextStyle`] naming any other family renders only after the caller
//! registers a face for it with `plotters::style::register_font`; otherwise
//! plotters reports the family as unavailable and the drawing call returns a
//! [`VisualizationError`].
//!
//! [`TextStyle`]: crate::visualizations::TextStyle

use crate::error::{VisualizationError, VisualizationResult};
use plotters::style::{FontStyle, register_font};
use std::sync::OnceLock;

/// Family name under which the embedded face is registered.
pub(super) const EMBEDDED_FONT_FAMILY: &str = "sans-serif";

static REGULAR_FACE: &[u8] = include_bytes!("../../../fonts/DejaVuSans.ttf");

/// Outcome of the one process-wide registration; `false` records a face that
/// failed to parse.
static REGISTRATION: OnceLock<bool> = OnceLock::new();

/// Register the embedded face under [`EMBEDDED_FONT_FAMILY`] once per process.
///
/// Every text-drawing entry point of this backend calls this before laying
/// out text, so no caller-visible initialization order exists.
///
/// # Errors
///
/// Returns [`VisualizationError::RenderingError`] when the embedded face
/// fails to parse as an OpenType font.
pub(super) fn register_embedded_typeface() -> VisualizationResult<()> {
    let registered = *REGISTRATION.get_or_init(|| {
        register_font(EMBEDDED_FONT_FAMILY, FontStyle::Normal, REGULAR_FACE).is_ok()
    });
    if registered {
        Ok(())
    } else {
        Err(VisualizationError::rendering_error(
            "embedded DejaVu Sans face failed to parse as an OpenType font",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{EMBEDDED_FONT_FAMILY, register_embedded_typeface};
    use plotters::style::{FontDesc, FontError, FontStyle};

    #[test]
    fn embedded_face_advance_scales_linearly_with_size() {
        register_embedded_typeface().expect("embedded face must parse");
        let width_at = |size: f64| {
            FontDesc::new(EMBEDDED_FONT_FAMILY.into(), size, FontStyle::Normal)
                .box_size("Hyperlegible")
                .expect("registered family must lay out")
                .0
        };
        let narrow = width_at(20.0);
        let wide = width_at(40.0);
        // Glyph advances scale linearly with pixel size and plotters truncates
        // the summed advance once, so floor(2x) lies in {2 floor(x), 2 floor(x) + 1}.
        assert!(narrow > 0, "a visible string must advance");
        assert!(
            wide == 2 * narrow || wide == 2 * narrow + 1,
            "40 px width {wide} must be twice the 20 px width {narrow} up to truncation"
        );
    }

    #[test]
    fn unregistered_family_is_reported_unavailable() {
        register_embedded_typeface().expect("embedded face must parse");
        let font = FontDesc::new(
            "cfd-schematics-unregistered".into(),
            20.0,
            FontStyle::Normal,
        );
        assert!(matches!(
            font.box_size("text"),
            Err(FontError::FontUnavailable)
        ));
    }
}
