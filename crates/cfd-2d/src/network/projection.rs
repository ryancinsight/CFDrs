//! Schematics-driven projection of routed network channels into 2D solver grids.
//!
//! The projection layer converts authoritative `cfd-schematics` paths and
//! cross-section metadata into fluid masks, solver-domain extents, and
//! per-channel occupancy summaries. This keeps the 2D solver aligned with the
//! topology authority instead of assuming rectangular channels.

use cfd_core::error::{Error, Result as CfdResult};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use cfd_schematics::domain::model::{ChannelSpec, CrossSectionSpec};
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;

use crate::solvers::ns_fvm::NavierStokesSolver2D;

use super::types::ChannelProjectionSummary;

#[derive(Debug, Clone, Copy)]
pub(crate) struct ProjectionDomain {
    pub(crate) length_m: f64,
    pub(crate) width_m: f64,
}

#[derive(Debug, Clone, Copy)]
struct PathMetrics {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    length: f64,
}

/// Aggregate projection summary for a full blueprint.
///
/// The summary is computed after every channel has been rasterized into the
/// 2D solver grid and captures both occupancy and per-channel projection
/// metadata.
#[derive(Debug, Clone)]
pub struct NetworkProjectionSummary<T> {
    /// Per-channel projection summaries in blueprint order.
    pub channel_summaries: Vec<ChannelProjectionSummary<T>>,
    /// Total number of fluid cells flagged across all channels.
    pub total_fluid_cell_count: usize,
    /// Mean fluid occupancy fraction across all projected channels.
    pub mean_fluid_fraction: T,
}

impl<T> NetworkProjectionSummary<T> {
    /// Number of projected channels in the summary.
    #[must_use]
    pub fn channel_count(&self) -> usize {
        self.channel_summaries.len()
    }
}

fn path_metrics(path: &[(f64, f64)]) -> Option<PathMetrics> {
    let first = path.first()?;
    let mut x_min = first.0;
    let mut x_max = first.0;
    let mut y_min = first.1;
    let mut y_max = first.1;
    let mut length = 0.0_f64;

    for window in path.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        x_min = x_min.min(x0).min(x1);
        x_max = x_max.max(x0).max(x1);
        y_min = y_min.min(y0).min(y1);
        y_max = y_max.max(y0).max(y1);
        length += (x1 - x0).hypot(y1 - y0);
    }

    Some(PathMetrics {
        x_min,
        x_max,
        y_min,
        y_max,
        length,
    })
}

fn venturi_width_at_x(venturi: &VenturiGeometryMetadata, x: f64, total_length: f64) -> f64 {
    let inlet = venturi.inlet_width_m.max(1e-18);
    let throat = venturi.throat_width_m.max(1e-18);
    let outlet = venturi.outlet_width_m.max(1e-18);
    let l_throat = venturi.throat_length_m.max(0.0);
    let remaining = (total_length - l_throat).max(0.0);
    let l_converge = remaining * venturi.throat_position.clamp(0.0, 1.0);
    let l_diverge = (remaining - l_converge).max(0.0);

    let x_inlet_end = 0.0;
    let x_converge_end = l_converge;
    let x_throat_end = x_converge_end + l_throat;
    let x_diverge_end = x_throat_end + l_diverge;

    if x <= x_inlet_end || l_converge <= 0.0 {
        inlet
    } else if x < x_converge_end {
        let frac = ((x - x_inlet_end) / l_converge).clamp(0.0, 1.0);
        inlet + frac * (throat - inlet)
    } else if x < x_throat_end {
        throat
    } else if x < x_diverge_end && l_diverge > 0.0 {
        let frac = ((x - x_throat_end) / l_diverge).clamp(0.0, 1.0);
        throat + frac * (outlet - throat)
    } else {
        outlet
    }
}

fn interpolate_piecewise_linear(points: &[(f64, f64)], x: f64) -> f64 {
    if points.is_empty() {
        return 0.0;
    }
    if points.len() == 1 || x <= points[0].0 {
        return points[0].1;
    }
    if x >= points[points.len() - 1].0 {
        return points[points.len() - 1].1;
    }

    for window in points.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        if x >= x0 && x <= x1 {
            let span = (x1 - x0).max(1e-18);
            let frac = ((x - x0) / span).clamp(0.0, 1.0);
            return y0 + frac * (y1 - y0);
        }
    }

    points[points.len() - 1].1
}

fn transformed_centerline(
    channel_id: &str,
    path: &[(f64, f64)],
    metrics: PathMetrics,
    length_scale: f64,
    width_scale: f64,
    y_offset: f64,
) -> CfdResult<Vec<(f64, f64)>> {
    if path.len() < 2 || metrics.length <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "Channel '{channel_id}' has a degenerate routed path for 2D projection",
        )));
    }

    let mut transformed = Vec::with_capacity(path.len());
    let mut arclength = 0.0_f64;
    transformed.push((0.0, y_offset + (path[0].1 - metrics.y_min) * width_scale));

    for window in path.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        arclength += (x1 - x0).hypot(y1 - y0);
        transformed.push((
            arclength * length_scale,
            y_offset + (y1 - metrics.y_min) * width_scale,
        ));
    }

    Ok(transformed)
}

fn default_half_width(channel: &ChannelSpec) -> f64 {
    match channel.cross_section {
        CrossSectionSpec::Circular { diameter_m } => diameter_m * 0.5,
        CrossSectionSpec::Rectangular { width_m, .. } => width_m * 0.5,
    }
}

#[must_use]
pub(crate) fn channel_projection_domain(channel: &ChannelSpec) -> ProjectionDomain {
    let base_width = channel
        .venturi_geometry
        .as_ref()
        .map_or_else(
            || default_half_width(channel) * 2.0,
            |geom| {
                geom.inlet_width_m
                    .max(geom.throat_width_m)
                    .max(geom.outlet_width_m)
            },
        )
        .max(default_half_width(channel) * 2.0);

    let metrics = path_metrics(&channel.path).unwrap_or(PathMetrics {
        x_min: 0.0,
        x_max: channel.length_m,
        y_min: 0.0,
        y_max: base_width,
        length: channel.length_m,
    });

    let width_scale = if metrics.length > 0.0 {
        channel.length_m / metrics.length
    } else {
        1.0
    };
    let layout_span_y = (metrics.y_max - metrics.y_min).max(0.0) * width_scale;
    let pad = base_width.max(1e-6);
    let length_m = channel
        .length_m
        .max((metrics.x_max - metrics.x_min) * width_scale)
        .max(1e-6);
    let width_m = ((layout_span_y * 0.1) + 2.0 * pad).max(1e-6);

    ProjectionDomain { length_m, width_m }
}

pub(crate) fn populate_channel_projection_mask<T>(
    solver: &mut NavierStokesSolver2D<T>,
    channel: &ChannelSpec,
    domain: ProjectionDomain,
) -> CfdResult<ChannelProjectionSummary<T>>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let metrics = path_metrics(&channel.path).ok_or_else(|| {
        Error::InvalidInput(format!(
            "Channel '{}' path metrics unavailable for projection",
            channel.id.as_str()
        ))
    })?;
    let length_scale = if metrics.length > 0.0 {
        domain.length_m / metrics.length
    } else {
        1.0
    };

    let y_span = (metrics.y_max - metrics.y_min).max(0.0);
    let pad = default_half_width(channel).max(1e-6);
    let y_scale = if y_span > 0.0 {
        (domain.width_m - 2.0 * pad) / y_span
    } else {
        1.0
    };
    let y_offset = (domain.width_m - y_span * y_scale) * 0.5;

    let transformed = transformed_centerline(
        channel.id.as_str(),
        &channel.path,
        metrics,
        length_scale,
        y_scale,
        y_offset,
    )?;
    let fluid_half_width = default_half_width(channel).max(1e-18);
    let venturi_meta = channel.venturi_geometry.clone();

    let nx = solver.grid.nx;
    let ny = solver.grid.ny;
    let grid_half_step = if ny > 0 {
        (domain.width_m / ny as f64) * 0.5
    } else {
        domain.width_m * 0.5
    }
    .max(1e-18);
    let mut fluid_cells = 0usize;

    for i in 0..nx {
        let x = solver.grid.x_center(i).to_f64().unwrap_or(0.0);
        let center_y = interpolate_piecewise_linear(&transformed, x);
        let half_width = match venturi_meta.as_ref() {
            Some(geom) => venturi_width_at_x(geom, x, domain.length_m).max(grid_half_step),
            None => fluid_half_width.max(grid_half_step),
        };

        for j in 0..ny {
            let y = solver.grid.y_center(j).to_f64().unwrap_or(0.0);
            let is_fluid = (y - center_y).abs() <= half_width;
            solver.field.mask[(i, j)] = is_fluid;
            if is_fluid {
                fluid_cells += 1;
            }
        }
    }

    let total_cells = nx.saturating_mul(ny).max(1);
    let fluid_fraction = T::from_f64(fluid_cells as f64 / total_cells as f64)
        .expect("analytical constant conversion");

    Ok(ChannelProjectionSummary {
        channel_id: channel.id.as_str().to_owned(),
        grid_length_m: T::from_f64(domain.length_m).expect("analytical constant conversion"),
        grid_width_m: T::from_f64(domain.width_m).expect("analytical constant conversion"),
        path_length_m: T::from_f64(metrics.length).expect("analytical constant conversion"),
        path_span_x_m: T::from_f64(metrics.x_max - metrics.x_min)
            .expect("analytical constant conversion"),
        path_span_y_m: T::from_f64(y_span).expect("analytical constant conversion"),
        fluid_cell_count: fluid_cells,
        fluid_fraction,
    })
}

pub(crate) fn summarize_projection<T>(
    summaries: Vec<ChannelProjectionSummary<T>>,
) -> NetworkProjectionSummary<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    let total_fluid_cell_count = summaries
        .iter()
        .map(|summary| summary.fluid_cell_count)
        .sum();
    let mean_fluid_fraction = if summaries.is_empty() {
        T::zero()
    } else {
        let sum = summaries
            .iter()
            .fold(T::zero(), |acc, summary| acc + summary.fluid_fraction);
        sum / T::from_usize(summaries.len()).expect("analytical constant conversion")
    };

    NetworkProjectionSummary {
        channel_summaries: summaries,
        total_fluid_cell_count,
        mean_fluid_fraction,
    }
}
