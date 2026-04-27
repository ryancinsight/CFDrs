/// Compute the optimal box dimensions for a given split pattern.
#[must_use]
pub fn adaptive_box_dims(
    plate_width_mm: f64,
    plate_height_mm: f64,
    total_branches: usize,
    channel_width_mm: f64,
    wall_clearance_mm: f64,
) -> (f64, f64) {
    let n = total_branches.max(1) as f64;
    let min_height = n * channel_width_mm + (n + 1.0) * wall_clearance_mm;

    let padding_factor = if total_branches <= 4 {
        1.8
    } else if total_branches <= 16 {
        1.5
    } else {
        1.2
    };

    let desired_height = (min_height * padding_factor).max(plate_height_mm * 0.3);
    let height = desired_height.min(plate_height_mm);
    let width = plate_width_mm;

    (width, height)
}
