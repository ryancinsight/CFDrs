//! Color-ID encoding/decoding for GPU pick pass.
//!
//! Encodes `(node_index, face_index)` into 4 bytes (RGBA) and decodes them
//! back. The encoding uses 16 bits per index (supports up to 65535 nodes
//! and 65535 faces per node).

/// Encode a (node_index, sub_element_index) pair into RGBA bytes.
///
/// Layout: R = node[7:0], G = node[15:8], B = sub[7:0], A = sub[15:8].
#[must_use]
pub fn encode(node_idx: u32, sub_idx: u32) -> [u8; 4] {
    [
        (node_idx & 0xFF) as u8,
        ((node_idx >> 8) & 0xFF) as u8,
        (sub_idx & 0xFF) as u8,
        ((sub_idx >> 8) & 0xFF) as u8,
    ]
}

/// Decode RGBA bytes back to (node_index, sub_element_index).
///
/// Returns `None` for the background (all zeros).
#[must_use]
pub fn decode(rgba: [u8; 4]) -> Option<(u32, u32)> {
    let node = rgba[0] as u32 | ((rgba[1] as u32) << 8);
    let sub = rgba[2] as u32 | ((rgba[3] as u32) << 8);
    if node == 0 && sub == 0 {
        return None;
    }
    Some((node, sub))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_small_values() {
        let encoded = encode(5, 10);
        assert_eq!(decode(encoded), Some((5, 10)));
    }

    #[test]
    fn round_trip_large_values() {
        let encoded = encode(1000, 50000);
        assert_eq!(decode(encoded), Some((1000, 50000)));
    }

    #[test]
    fn zero_is_background() {
        assert_eq!(decode([0, 0, 0, 0]), None);
    }

    #[test]
    fn node_one_face_zero_is_valid() {
        let encoded = encode(1, 0);
        assert_eq!(decode(encoded), Some((1, 0)));
    }
}
