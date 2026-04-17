use std::collections::HashMap;

use cfd_core::error::Result as CfdResult;
use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{ChannelShape, NetworkBlueprint};
use cfd_schematics::domain::rules::BlueprintValidator;
use cfd_schematics::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};

use super::projection::{
    channel_projection_domain, populate_channel_projection_mask, summarize_projection,
};
use super::reference::solve_reference_trace;
use super::types::{Channel2dEntry, Network2DSolver};
use super::validate_blueprint_for_2d_projection;
use super::ChannelReferenceTrace;

/// A [`GraphSink`] that converts a validated [`NetworkBlueprint`] into a
/// solver-ready [`Network2DSolver<T>`].
pub struct Network2dBuilderSink<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    blood: BloodModel<T>,
    density: f64,
    /// Total inlet flow rate [m³/s] used to scale the cfd-1d reference trace.
    total_flow_rate_m3_s: f64,
    grid_nx: usize,
    grid_ny: usize,
    separation_tracking_enabled: bool,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> Network2dBuilderSink<T> {
    /// Create a new sink.
    #[must_use]
    pub fn new(
        blood: BloodModel<T>,
        density: f64,
        total_flow_rate_m3_s: f64,
        grid_nx: usize,
        grid_ny: usize,
    ) -> Self {
        Self {
            blood,
            density,
            total_flow_rate_m3_s,
            grid_nx,
            grid_ny,
            separation_tracking_enabled: false,
        }
    }

    /// Enable or disable the expensive cell-tracking separation postprocess.
    #[must_use]
    pub fn with_separation_tracking(mut self, enabled: bool) -> Self {
        self.separation_tracking_enabled = enabled;
        self
    }
}

impl<T> GraphSink for Network2dBuilderSink<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug + 'static,
{
    type Output = Network2DSolver<T>;

    fn build(&self, blueprint: &NetworkBlueprint) -> CfdResult<Network2DSolver<T>> {
        BlueprintValidator::validate(blueprint)?;
        validate_blueprint_for_2d_projection(blueprint)?;

        let mu = blood_viscosity_f64(&self.blood);
        let reference_trace = solve_reference_trace::<T>(
            blueprint,
            self.density,
            mu,
            self.total_flow_rate_m3_s,
        )?;
        let channel_reference_by_id: HashMap<&str, &ChannelReferenceTrace<T>> = reference_trace
            .channel_traces
            .iter()
            .map(|trace| (trace.channel_id.as_str(), trace))
            .collect();

        let mut entries = Vec::with_capacity(blueprint.channels.len());
        let mut projection_summaries = Vec::with_capacity(blueprint.channels.len());

        for channel in &blueprint.channels {
            let channel_reference = channel_reference_by_id
                .get(channel.id.as_str())
                .copied()
                .ok_or_else(|| {
                    cfd_core::error::Error::InvalidInput(format!(
                        "Network2DSolver missing cfd-1d reference trace for channel '{}'",
                        channel.id.as_str()
                    ))
                })?;
            let cross_section_area_m2 = channel.cross_section.area();
            let q_ch = channel_reference.flow_rate_m3_s.to_f64().unwrap_or(0.0);
            let mean_velocity_m_s = q_ch / cross_section_area_m2.max(1e-18);

            let therapy_zone = channel
                .therapy_zone
                .or_else(|| {
                    channel
                        .metadata
                        .as_ref()
                        .and_then(|m| m.get::<TherapyZoneMetadata>())
                        .map(|tz| tz.zone)
                })
                .unwrap_or(TherapyZone::MixedFlow);

            let venturi_meta = channel.venturi_geometry.clone().or_else(|| {
                channel
                    .metadata
                    .as_ref()
                    .and_then(|m| m.get::<VenturiGeometryMetadata>())
                    .cloned()
            });
            let is_venturi_throat = venturi_meta.is_some();

            let projection_domain = channel_projection_domain(channel);
            let l_t = T::from_f64(projection_domain.length_m).expect("analytical constant conversion");
            let grid_w = T::from_f64(projection_domain.width_m)
                .expect("analytical constant conversion");
            let grid = StaggeredGrid2D::new(self.grid_nx, self.grid_ny, l_t, grid_w);
            let density_t = T::from_f64(self.density).expect("analytical constant conversion");
            let config = solver_config_for_channel::<T>(&channel.channel_shape, mean_velocity_m_s);
            let mut solver = NavierStokesSolver2D::new(grid, self.blood.clone(), density_t, config);

            let projection =
                populate_channel_projection_mask(&mut solver, channel, projection_domain)?;
            projection_summaries.push(projection.clone());

            entries.push(Channel2dEntry {
                id: channel.id.as_str().to_owned(),
                therapy_zone,
                is_venturi_throat,
                projection,
                flow_rate_m3_s: q_ch,
                cross_section: channel.cross_section,
                cross_section_area_m2,
                length_m: projection_domain.length_m,
                viscosity_pa_s: mu,
                reference_trace: channel_reference.clone(),
                solver,
            });
        }

        let projection = summarize_projection(projection_summaries);

        Ok(Network2DSolver {
            blueprint: blueprint.clone(),
            channels: entries,
            reference_trace,
            projection,
            reference_density_kg_m3: self.density,
            reference_viscosity_pa_s: mu,
            separation_tracking_enabled: self.separation_tracking_enabled,
        })
    }
}

/// Extract a representative dynamic viscosity [Pa·s] from the blood model.
fn blood_viscosity_f64<T: RealField + Copy + Float + FromPrimitive + ToPrimitive>(
    model: &BloodModel<T>,
) -> f64 {
    const REF_SHEAR: f64 = 100.0;
    let shear_t = T::from_f64(REF_SHEAR).expect("analytical constant conversion");
    model.viscosity(shear_t).to_f64().unwrap_or(3.5e-3)
}

/// Select a SIMPLE/PISO relaxation profile for the channel shape.
fn solver_config_for_channel<T: RealField + Copy + Float + FromPrimitive>(
    shape: &ChannelShape,
    mean_velocity_m_s: f64,
) -> SIMPLEConfig<T> {
    if mean_velocity_m_s > 0.55 {
        SIMPLEConfig {
            max_iterations: 12_000,
            alpha_u: T::from_f64(0.25).expect("analytical constant conversion"),
            alpha_p: T::from_f64(0.05).expect("analytical constant conversion"),
            n_correctors: 4,
            ..SIMPLEConfig::default()
        }
    } else if matches!(shape, ChannelShape::Serpentine { .. }) || mean_velocity_m_s > 0.5 {
        SIMPLEConfig {
            max_iterations: 8_000,
            alpha_u: T::from_f64(0.45).expect("analytical constant conversion"),
            alpha_p: T::from_f64(0.15).expect("analytical constant conversion"),
            n_correctors: 2,
            ..SIMPLEConfig::default()
        }
    } else {
        SIMPLEConfig {
            max_iterations: 5_000,
            alpha_u: T::from_f64(0.5).expect("analytical constant conversion"),
            alpha_p: T::from_f64(0.2).expect("analytical constant conversion"),
            ..SIMPLEConfig::default()
        }
    }
}
