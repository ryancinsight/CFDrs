//! Mode toolbar — select, move, rotate, scale modes.

/// The active interaction mode in the viewport.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum EditMode {
    /// Object selection mode.
    #[default]
    Select,
    /// Object translation mode.
    Move,
    /// Object rotation mode.
    Rotate,
    /// Object scaling mode.
    Scale,
    /// Measurement mode (pick geometry to measure).
    Measure,
    /// 2D sketch drawing mode.
    Sketch,
}
