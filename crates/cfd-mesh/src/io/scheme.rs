//! Scheme import: read millifluidic chip designs from cfd-schematics JSON.
//!
//! Provides two import paths:
//!
//! 1. **Raw JSON** — `import_schematic()` parses a simple JSON format with
//!    `substrate` and `channels` keys (no dependency on `cfd-schematics`).
//!
//! 2. **`cfd-schematics` interchange** — `from_channel_system()` converts
//!    a `cfd_schematics::ChannelSystem` into a `Schematic` that the
//!    channel/sweep pipeline can mesh directly.

use std::io::Read;

use crate::core::scalar::{Real, Point3r};
use crate::core::error::{MeshError, MeshResult};
use crate::channel::profile::ChannelProfile;
use crate::channel::path::ChannelPath;

/// A parsed channel definition from a schematic.
#[derive(Clone, Debug)]
pub struct ChannelDef {
    /// Channel identifier.
    pub id: String,
    /// Channel path (centerline).
    pub path: ChannelPath,
    /// Channel cross-section profile.
    pub profile: ChannelProfile,
    /// Optional width scaling factors along the path (for variable-width channels).
    pub width_scales: Option<Vec<Real>>,
}

/// A parsed substrate definition.
#[derive(Clone, Debug)]
pub struct SubstrateDef {
    /// Width (X).
    pub width: Real,
    /// Depth (Y).
    pub depth: Real,
    /// Height (Z).
    pub height: Real,
    /// Origin.
    pub origin: Point3r,
}

/// A fully parsed millifluidic schematic.
#[derive(Clone, Debug)]
pub struct Schematic {
    /// Substrate definition.
    pub substrate: SubstrateDef,
    /// Channel definitions.
    pub channels: Vec<ChannelDef>,
}

/// Import a schematic from a JSON reader.
///
/// Expected JSON format:
/// ```json
/// {
///   "substrate": { "width": 30.0, "depth": 20.0, "height": 10.0 },
///   "channels": [
///     {
///       "id": "ch1",
///       "diameter": 1.0,
///       "segments": 16,
///       "path": [[0,0,5], [10,0,5], [20,0,5]]
///     }
///   ]
/// }
/// ```
pub fn import_schematic<R: Read>(reader: R) -> MeshResult<Schematic> {
    let value: serde_json::Value =
        serde_json::from_reader(reader).map_err(MeshError::Json)?;

    let substrate = parse_substrate(&value)?;
    let channels = parse_channels(&value)?;

    Ok(Schematic {
        substrate,
        channels,
    })
}

fn parse_substrate(value: &serde_json::Value) -> MeshResult<SubstrateDef> {
    let sub = value
        .get("substrate")
        .ok_or_else(|| MeshError::Other("missing 'substrate' key".to_string()))?;

    let width = sub["width"]
        .as_f64()
        .ok_or_else(|| MeshError::Other("missing substrate width".to_string()))?
        as Real;
    let depth = sub["depth"]
        .as_f64()
        .ok_or_else(|| MeshError::Other("missing substrate depth".to_string()))?
        as Real;
    let height = sub["height"]
        .as_f64()
        .ok_or_else(|| MeshError::Other("missing substrate height".to_string()))?
        as Real;

    let origin = if let Some(o) = sub.get("origin") {
        Point3r::new(
            o[0].as_f64().unwrap_or(0.0) as Real,
            o[1].as_f64().unwrap_or(0.0) as Real,
            o[2].as_f64().unwrap_or(0.0) as Real,
        )
    } else {
        Point3r::origin()
    };

    Ok(SubstrateDef {
        width,
        depth,
        height,
        origin,
    })
}

fn parse_channels(value: &serde_json::Value) -> MeshResult<Vec<ChannelDef>> {
    let channels = value
        .get("channels")
        .and_then(|c| c.as_array())
        .ok_or_else(|| MeshError::Other("missing 'channels' array".to_string()))?;

    let mut defs = Vec::with_capacity(channels.len());

    for ch in channels {
        let id = ch["id"]
            .as_str()
            .unwrap_or("unnamed")
            .to_string();

        let diameter = ch["diameter"]
            .as_f64()
            .ok_or_else(|| MeshError::Other("missing channel diameter".to_string()))?
            as Real;
        let segments = ch["segments"].as_u64().unwrap_or(16) as usize;

        let path_arr = ch["path"]
            .as_array()
            .ok_or_else(|| MeshError::Other("missing channel path".to_string()))?;

        let points: Vec<Point3r> = path_arr
            .iter()
            .map(|p| {
                let arr = p.as_array().ok_or_else(|| {
                    MeshError::Other("path point must be array".to_string())
                })?;
                Ok(Point3r::new(
                    arr[0].as_f64().unwrap_or(0.0) as Real,
                    arr[1].as_f64().unwrap_or(0.0) as Real,
                    arr[2].as_f64().unwrap_or(0.0) as Real,
                ))
            })
            .collect::<MeshResult<Vec<_>>>()?;

        defs.push(ChannelDef {
            id,
            path: ChannelPath::new(points),
            profile: ChannelProfile::Circular {
                radius: diameter / 2.0,
                segments,
            },
            width_scales: None,
        });
    }

    Ok(defs)
}

// =========================================================================
// cfd-schematics bridge (behind `scheme-io` feature)
// =========================================================================

/// Convert a `cfd_schematics::ChannelSystem` into a `Schematic`.
///
/// 2D channel centerlines are lifted to 3D at `z = height / 2` so channels
/// run through the centre of the substrate.
///
/// # Arguments
/// - `system` — the channel system from `cfd_schematics::geometry::generator::create_geometry`
/// - `height` — substrate height in mm (also used as z-centering)
/// - `channel_segments` — number of cross-section segments per channel
pub fn from_channel_system(
    system: &cfd_schematics::geometry::ChannelSystem,
    height: Real,
    channel_segments: usize,
) -> MeshResult<Schematic> {
    let interchange = system.to_interchange();
    from_interchange(&interchange, height, channel_segments)
}

/// Convert an `InterchangeChannelSystem` into a `Schematic`.
pub fn from_interchange(
    interchange: &cfd_schematics::geometry::InterchangeChannelSystem,
    height: Real,
    channel_segments: usize,
) -> MeshResult<Schematic> {
    let (bw, bd) = interchange.box_dims_mm;

    let substrate = SubstrateDef {
        width: bw as Real,
        depth: bd as Real,
        height,
        origin: Point3r::origin(),
    };

    let mut channels = Vec::with_capacity(interchange.channels.len());

    for ch in &interchange.channels {
        if ch.centerline_mm.len() < 2 {
            continue;
        }

        let mid_z = height / 2.0;

        // Lift 2D centerline to 3D at z = mid_z
        let points: Vec<Point3r> = ch
            .centerline_mm
            .iter()
            .map(|&(x, y)| Point3r::new(x as Real, y as Real, mid_z))
            .collect();

        let mut width_scales = None;

        let profile = match &ch.profile {
            cfd_schematics::geometry::InterchangeChannelProfile::Constant {
                width_mm,
                height_mm,
                ..
            } => ChannelProfile::Rectangular {
                width: *width_mm as Real,
                height: *height_mm as Real,
            },
            cfd_schematics::geometry::InterchangeChannelProfile::Frustum {
                inlet_width_mm,
                height_mm,
                width_profile_mm,
                ..
            } => {
                let inlet_w = *inlet_width_mm as Real;
                // Calculate scaling factors relative to inlet width
                let scales: Vec<Real> = width_profile_mm
                    .iter()
                    .map(|&w| (w as Real) / inlet_w)
                    .collect();
                
                width_scales = Some(scales);

                ChannelProfile::Rectangular {
                    width: inlet_w,
                    height: *height_mm as Real,
                }
            }
        };

        channels.push(ChannelDef {
            id: format!("ch{}", ch.id),
            path: ChannelPath::new(points),
            profile,
            width_scales,
        });
    }

    if channels.is_empty() {
        return Err(MeshError::ChannelError {
            message: "no channels with >= 2 centerline points".to_string(),
        });
    }

    Ok(Schematic {
        substrate,
        channels,
    })
}
