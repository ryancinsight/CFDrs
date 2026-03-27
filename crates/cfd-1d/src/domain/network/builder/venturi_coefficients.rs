//! Venturi-specific resistance coefficient computation.
//!
//! Extracts venturi tube coefficients from schematics geometry metadata for
//! use in the network builder pipeline.

use cfd_core::error::Result;
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Compute venturi (R, K) pair from schematics geometry metadata.
///
/// Uses the analytical Venturi model with discharge and expansion coefficients
/// derived from the convergent/divergent half-angles.
pub(crate) fn venturi_coefficients<T, F>(
    metadata: &cfd_schematics::geometry::metadata::VenturiGeometryMetadata,
    fluid: &F,
) -> Result<(T, T)>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    use crate::physics::resistance::models::{
        ExpansionType, FlowConditions as ModelFlowConditions, ResistanceModel, VenturiGeometry,
        VenturiModel,
    };
    use crate::{discharge_coefficient_from_convergent_half_angle_deg, venturi_taper_length_m};

    let hydraulic_diameter =
        |width_m: f64, height_m: f64| 2.0 * width_m * height_m / (width_m + height_m).max(1e-18);
    let inlet_d = hydraulic_diameter(metadata.inlet_width_m, metadata.throat_height_m);
    let throat_d = hydraulic_diameter(metadata.throat_width_m, metadata.throat_height_m);
    let outlet_d = hydraulic_diameter(metadata.outlet_width_m, metadata.throat_height_m);
    let conv_len = venturi_taper_length_m(
        metadata.inlet_width_m,
        metadata.throat_width_m,
        metadata.convergent_half_angle_deg,
    );
    let diff_len = venturi_taper_length_m(
        metadata.outlet_width_m,
        metadata.throat_width_m,
        metadata.divergent_half_angle_deg,
    );
    let total_length = metadata.throat_length_m + conv_len + diff_len;

    let inlet_d_t = T::from_f64(inlet_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration("inlet diameter conversion failed".to_string())
    })?;
    let throat_d_t = T::from_f64(throat_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "throat diameter conversion failed".to_string(),
        )
    })?;
    let outlet_d_t = T::from_f64(outlet_d).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "outlet diameter conversion failed".to_string(),
        )
    })?;
    let throat_len_t = T::from_f64(total_length).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration("throat length conversion failed".to_string())
    })?;
    let total_len_t = T::from_f64(total_length).ok_or_else(|| {
        cfd_core::error::Error::InvalidConfiguration(
            "venturi total length conversion failed".to_string(),
        )
    })?;

    let mut model = VenturiModel::new(inlet_d_t, throat_d_t, outlet_d_t, throat_len_t, total_len_t)
        .with_geometry(VenturiGeometry::Custom {
            discharge_coefficient: discharge_coefficient_from_convergent_half_angle_deg(
                metadata.convergent_half_angle_deg,
            ),
        })
        .with_expansion(ExpansionType::Gradual {
            half_angle_deg: metadata.divergent_half_angle_deg.clamp(1.0, 45.0),
        });
    model.throat_roughness =
        T::from_f64(1e-7).expect("Mathematical constant conversion compromised");

    let mut conds = ModelFlowConditions::new(T::zero());
    conds.flow_rate =
        Some(T::from_f64(1e-12).expect("Mathematical constant conversion compromised"));
    model.calculate_coefficients(fluid, &conds)
}
