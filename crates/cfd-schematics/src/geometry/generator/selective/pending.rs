use super::super::super::types::Point2D;

#[derive(Debug, Clone, Copy)]
pub(super) struct PendingVenturiPath {
    pub(super) channel_idx: usize,
    pub(super) start: Option<Point2D>,
    pub(super) end: Option<Point2D>,
    pub(super) preferred_y: f64,
    pub(super) fallback_length_m: f64,
}

