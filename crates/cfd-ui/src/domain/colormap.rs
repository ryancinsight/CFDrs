//! Region and field colormaps for mesh visualization.

/// A colormap for mapping scalar values to RGBA colors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Colormap {
    /// Blue-to-red diverging colormap.
    BlueRed,
    /// Viridis perceptually uniform colormap.
    Viridis,
    /// Grayscale.
    Grayscale,
    /// Distinct colors per region (qualitative palette).
    RegionPalette,
}

/// An RGBA color with 8-bit channels.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Color {
    /// Create an opaque color.
    #[must_use]
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: 255 }
    }

    /// Convert to `[f32; 4]` normalized to [0, 1].
    #[must_use]
    pub fn to_f32_array(self) -> [f32; 4] {
        [
            self.r as f32 / 255.0,
            self.g as f32 / 255.0,
            self.b as f32 / 255.0,
            self.a as f32 / 255.0,
        ]
    }
}

/// Map a normalized value `t` in [0, 1] to an RGBA color using the given colormap.
///
/// # Theorem — Continuity
///
/// For the `BlueRed` and `Viridis` colormaps, the output color is a
/// piecewise-linear function of `t`, guaranteeing C^0 continuity across
/// the entire [0, 1] domain.  **Proof sketch**: each channel is computed
/// as `lerp(c_low, c_high, t)` within a single segment, which is linear
/// in `t`.  ∎
#[must_use]
pub fn colorize(t: f64, colormap: Colormap) -> Color {
    let t = t.clamp(0.0, 1.0);
    match colormap {
        Colormap::BlueRed => {
            let r = (t * 255.0) as u8;
            let b = ((1.0 - t) * 255.0) as u8;
            let g = ((1.0 - (2.0 * t - 1.0).abs()) * 255.0) as u8;
            Color::rgb(r, g, b)
        }
        Colormap::Viridis => {
            // Simplified 3-stop Viridis approximation.
            let (r, g, b) = if t < 0.5 {
                let s = t * 2.0;
                (
                    lerp(68.0, 33.0, s),
                    lerp(1.0, 145.0, s),
                    lerp(84.0, 140.0, s),
                )
            } else {
                let s = (t - 0.5) * 2.0;
                (
                    lerp(33.0, 253.0, s),
                    lerp(145.0, 231.0, s),
                    lerp(140.0, 37.0, s),
                )
            };
            Color::rgb(r as u8, g as u8, b as u8)
        }
        Colormap::Grayscale => {
            let v = (t * 255.0) as u8;
            Color::rgb(v, v, v)
        }
        Colormap::RegionPalette => region_color(t),
    }
}

/// Map a region index (as normalized t = index / max) to a distinct color.
#[must_use]
pub fn region_color(t: f64) -> Color {
    let idx = (t * (REGION_PALETTE.len() - 1) as f64).round() as usize;
    let idx = idx.min(REGION_PALETTE.len() - 1);
    REGION_PALETTE[idx]
}

/// Linear interpolation.
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Qualitative palette for up to 10 distinct regions.
const REGION_PALETTE: [Color; 10] = [
    Color::rgb(31, 119, 180),   // blue
    Color::rgb(255, 127, 14),   // orange
    Color::rgb(44, 160, 44),    // green
    Color::rgb(214, 39, 40),    // red
    Color::rgb(148, 103, 189),  // purple
    Color::rgb(140, 86, 75),    // brown
    Color::rgb(227, 119, 194),  // pink
    Color::rgb(127, 127, 127),  // gray
    Color::rgb(188, 189, 34),   // olive
    Color::rgb(23, 190, 207),   // cyan
];
