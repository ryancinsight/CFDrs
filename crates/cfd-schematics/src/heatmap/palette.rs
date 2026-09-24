/// Interpolate hex colour: yellow (#FFD700) at 0.0, red (#CC0000) at 1.0.
pub(super) fn cancer_cav_color(v: f64) -> String {
    let t = v.clamp(0.0, 1.0);
    let r = ((0.8 + 0.2 * (1.0 - t)) * 255.0) as u8;
    let g = ((1.0 - t) * 215.0) as u8;
    format!("#{r:02X}{g:02X}00")
}

