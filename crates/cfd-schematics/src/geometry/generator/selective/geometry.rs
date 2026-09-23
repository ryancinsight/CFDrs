#[derive(Debug, Clone, Copy)]
pub(super) struct SelectiveTreeGeometry {
    pub(super) box_dims_mm: (f64, f64),
    pub(super) trunk_length_m: f64,
    pub(super) branch_length_m: f64,
    pub(super) hybrid_branch_length_m: f64,
    pub(super) main_width_m: f64,
    pub(super) throat_width_m: f64,
    pub(super) throat_length_m: f64,
    pub(super) channel_height_m: f64,
}

