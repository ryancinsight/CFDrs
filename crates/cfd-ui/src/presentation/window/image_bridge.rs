//! Converts wgpu offscreen BGRA pixel buffers to gpui `RenderImage`.
//!
//! # Theorem — Channel Swap Involution
//!
//! The BGRA-to-RGBA byte swap is its own inverse: applying it twice yields the
//! original buffer. **Proof sketch**: swapping indices 0 and 2 within each
//! 4-byte group is a transposition, and any transposition composed with itself
//! is the identity.  QED

use std::sync::Arc;
use gpui::RenderImage;
use smallvec::smallvec;

/// Convert a BGRA pixel buffer (from wgpu offscreen render) to an
/// `Arc<RenderImage>` suitable for display in a gpui `img()` element.
///
/// # Panics
///
/// Panics if `bgra.len() != (width * height * 4) as usize`.
pub fn bgra_to_render_image(bgra: Vec<u8>, width: u32, height: u32) -> Arc<RenderImage> {
    debug_assert_eq!(
        bgra.len(),
        (width as usize) * (height as usize) * 4,
        "pixel buffer length must equal width * height * 4",
    );

    let mut rgba = bgra;
    for pixel in rgba.chunks_exact_mut(4) {
        pixel.swap(0, 2);
    }

    let rgba_image = image::RgbaImage::from_raw(width, height, rgba)
        .expect("pixel buffer length equals width * height * 4");
    Arc::new(RenderImage::new(smallvec![image::Frame::new(rgba_image)]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swap_involution() {
        let original: Vec<u8> = vec![10, 20, 30, 255, 50, 60, 70, 128];
        let width = 2;
        let height = 1;

        // First swap: BGRA -> RGBA
        let mut data = original.clone();
        for pixel in data.chunks_exact_mut(4) {
            pixel.swap(0, 2);
        }
        // Second swap: RGBA -> BGRA (recovers original)
        for pixel in data.chunks_exact_mut(4) {
            pixel.swap(0, 2);
        }
        assert_eq!(data, original);

        // Ensure the function produces a valid image.
        let img = bgra_to_render_image(original, width, height);
        assert_eq!(img.frame_count(), 1);
    }
}
