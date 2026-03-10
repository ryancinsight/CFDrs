//! Scheme import: read millifluidic chip designs from cfd-schematics JSON.
//!
//! Provides two import paths:
//!
//! 1. **Raw JSON** — `import_schematic()` parses a simple JSON format with
//!    `substrate` and `channels` keys (no dependency on `cfd-schematics`).
//!
//! 2. **`cfd-schematics` interchange** — `from_blueprint()` converts
//!    a `cfd_schematics::domain::model::NetworkBlueprint` into a `Schematic` that the
//!    channel/sweep pipeline can mesh directly.

use std::io::Read;

use crate::application::channel::path::ChannelPath;
use crate::application::channel::profile::ChannelProfile;
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::scalar::{Point3r, Real};

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
    let value: serde_json::Value = serde_json::from_reader(reader).map_err(MeshError::Json)?;

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
        let id = ch["id"].as_str().unwrap_or("unnamed").to_string();

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
                let arr = p
                    .as_array()
                    .ok_or_else(|| MeshError::Other("path point must be array".to_string()))?;
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
// cfd-schematics bridge
// =========================================================================

/// Convert a `cfd_schematics::domain::model::NetworkBlueprint` into a `Schematic`.
///
/// 2D channel centerlines are lifted to 3D at `z = height / 2` so channels
/// run through the centre of the substrate.
///
/// # Arguments
/// - `blueprint` — the network blueprint from `cfd_schematics`
/// - `height` — substrate height in mm (also used as z-centering)
/// - `channel_segments` — number of cross-section segments per channel
pub fn from_blueprint(
    blueprint: &cfd_schematics::domain::model::NetworkBlueprint,
    height: Real,
    channel_segments: usize,
) -> MeshResult<Schematic> {
    let (bw, bd) = blueprint.box_dims;

    let substrate = SubstrateDef {
        width: bw as Real,
        depth: bd as Real,
        height,
        origin: Point3r::origin(),
    };

    let mut channels = Vec::with_capacity(blueprint.channels.len());
    let node_points: std::collections::HashMap<&str, (f64, f64)> = blueprint
        .nodes
        .iter()
        .map(|node| (node.id.as_str(), node.point))
        .collect();

    for ch in &blueprint.channels {
        let centerline_2d: Vec<(f64, f64)> = match ch.path.as_slice() {
            [] => node_points
                .get(ch.from.as_str())
                .zip(node_points.get(ch.to.as_str()))
                .map_or_else(Vec::new, |(&start, &end)| vec![start, end]),
            [midpoint] => {
                let mut points = Vec::with_capacity(3);
                if let Some(&start) = node_points.get(ch.from.as_str()) {
                    points.push(start);
                }
                points.push(*midpoint);
                if let Some(&end) = node_points.get(ch.to.as_str()) {
                    points.push(end);
                }
                points
            }
            points => points.to_vec(),
        };

        if centerline_2d.len() < 2 {
            continue;
        }

        let mid_z = height / 2.0;

        // Lift 2D centerline to 3D at z = mid_z
        let points: Vec<Point3r> = centerline_2d
            .iter()
            .map(|&(x, y)| Point3r::new(x as Real, y as Real, mid_z))
            .collect();

        let mut width_scales = None;

        let profile = if let Some(vg) = &ch.venturi_geometry {
            let inlet_w = (vg.inlet_width_m * 1000.0) as Real;
            let throat_w = (vg.throat_width_m * 1000.0) as Real;
            let outlet_w = (vg.outlet_width_m * 1000.0) as Real;

            // Reconstruct frustum width profile matching the lifted centerline.
            let n = centerline_2d.len();
            let mut scales = Vec::with_capacity(n);
            for idx in 0..n {
                let t = idx as f64 / (n - 1).max(1) as f64;
                let w = if t <= 0.5 {
                    let local = t / 0.5;
                    inlet_w + (throat_w - inlet_w) * local as Real
                } else {
                    let local = (t - 0.5) / 0.5;
                    throat_w + (outlet_w - throat_w) * local as Real
                };
                scales.push(w / inlet_w);
            }
            width_scales = Some(scales);

            let h = match ch.cross_section {
                cfd_schematics::domain::model::CrossSectionSpec::Rectangular {
                    height_m, ..
                } => (height_m * 1000.0) as Real,
                _ => inlet_w,
            };

            ChannelProfile::Rectangular {
                width: inlet_w,
                height: h,
            }
        } else {
            match ch.cross_section {
                cfd_schematics::domain::model::CrossSectionSpec::Circular { diameter_m } => {
                    ChannelProfile::Circular {
                        radius: (diameter_m as Real * 1000.0) / 2.0,
                        segments: channel_segments,
                    }
                }
                cfd_schematics::domain::model::CrossSectionSpec::Rectangular {
                    width_m,
                    height_m,
                } => ChannelProfile::Rectangular {
                    width: width_m as Real * 1000.0,
                    height: height_m as Real * 1000.0,
                },
            }
        };

        channels.push(ChannelDef {
            id: ch.id.as_str().to_string(),
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
