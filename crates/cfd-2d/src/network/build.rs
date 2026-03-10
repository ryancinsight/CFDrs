use std::collections::HashMap;

use cfd_core::error::Result as CfdResult;
use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{CrossSectionSpec, NetworkBlueprint};
use cfd_schematics::domain::rules::BlueprintValidator;
use cfd_schematics::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use crate::solvers::venturi_flow::VenturiGeometry;

use super::postprocess::{populate_circular_mask, populate_venturi_mask};
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
        }
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

        let reference_trace = solve_reference_trace::<T>(
            blueprint,
            self.density,
            blood_viscosity_f64(&self.blood),
            self.total_flow_rate_m3_s,
        )?;
        let channel_reference_by_id: HashMap<&str, &ChannelReferenceTrace<T>> = reference_trace
            .channel_traces
            .iter()
            .map(|trace| (trace.channel_id.as_str(), trace))
            .collect();
        let mu = blood_viscosity_f64(&self.blood);

        let mut entries = Vec::with_capacity(blueprint.channels.len());

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
            let q_ch = channel_reference.flow_rate_m3_s.to_f64().unwrap_or(0.0);

            let (w, _) = channel.cross_section.dims();
            let cross_section_area_m2 = channel.cross_section.area();

            let therapy_zone = channel
                .therapy_zone
                .clone()
                .or_else(|| {
                    channel
                        .metadata
                        .as_ref()
                        .and_then(|m| m.get::<TherapyZoneMetadata>())
                        .map(|tz| tz.zone.clone())
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

            let l = channel.length_m;
            let l_t = T::from_f64(l).unwrap_or_else(T::one);
            let grid_w = if let Some(vm) = venturi_meta.as_ref() {
                T::from_f64(vm.inlet_width_m).unwrap_or_else(T::one)
            } else {
                T::from_f64(w).unwrap_or_else(T::one)
            };
            let grid = StaggeredGrid2D::new(self.grid_nx, self.grid_ny, l_t, grid_w);
            let density_t = T::from_f64(self.density).unwrap_or_else(T::one);
            let mut config = SIMPLEConfig::default();
            config.max_iterations = 5_000;
            let mut solver = NavierStokesSolver2D::new(grid, self.blood.clone(), density_t, config);

            if let Some(vm) = venturi_meta.as_ref() {
                let converge_l = (l - vm.throat_length_m).max(0.0) / 2.0;
                let geom = VenturiGeometry::new(
                    T::from_f64(vm.inlet_width_m).unwrap_or_else(T::one),
                    T::from_f64(vm.throat_width_m).unwrap_or_else(T::one),
                    T::zero(),
                    T::from_f64(converge_l).unwrap_or_else(T::one),
                    T::from_f64(vm.throat_length_m).unwrap_or_else(T::one),
                    T::from_f64(converge_l).unwrap_or_else(T::one),
                    T::from_f64(vm.throat_height_m).unwrap_or_else(T::one),
                );
                populate_venturi_mask(&mut solver, &geom, self.grid_nx, self.grid_ny);
            } else if let CrossSectionSpec::Circular { diameter_m } = channel.cross_section {
                populate_circular_mask(&mut solver, diameter_m, self.grid_nx, self.grid_ny);
            }

            entries.push(Channel2dEntry {
                id: channel.id.as_str().to_owned(),
                therapy_zone,
                is_venturi_throat,
                flow_rate_m3_s: q_ch,
                cross_section: channel.cross_section,
                cross_section_area_m2,
                length_m: l,
                viscosity_pa_s: mu,
                reference_trace: channel_reference.clone(),
                solver,
            });
        }

        Ok(Network2DSolver {
            channels: entries,
            reference_trace,
        })
    }
}

/// Extract a representative dynamic viscosity [Pa·s] from the blood model.
fn blood_viscosity_f64<T: RealField + Copy + Float + FromPrimitive + ToPrimitive>(
    model: &BloodModel<T>,
) -> f64 {
    const REF_SHEAR: f64 = 100.0;
    let shear_t = T::from_f64(REF_SHEAR).unwrap_or_else(T::one);
    model.viscosity(shear_t).to_f64().unwrap_or(3.5e-3)
}
