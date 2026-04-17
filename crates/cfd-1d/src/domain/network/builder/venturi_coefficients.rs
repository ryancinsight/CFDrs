//! Venturi-specific resistance coefficient computation.
//!
//! Extracts venturi tube coefficients from schematics geometry metadata for
//! use in the network builder pipeline.
//!
//! The coefficient probe must be driven by an infinitesimal flow rate with the
//! explicit velocity left unset. Passing `velocity = Some(0)` forces the venturi
//! model into its zero-flow fallback and suppresses the angle-dependent inertial
//! loss coefficient.

use cfd_core::error::{Error, Result};
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

    let validate_positive = |name: &str, value: f64| -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(Error::InvalidConfiguration(format!(
                "Venturi {name} must be finite and positive"
            )))
        }
    };

    let validate_angle = |name: &str, value: f64| -> Result<f64> {
        if value.is_finite() && (1.0..=45.0).contains(&value) {
            Ok(value)
        } else {
            Err(Error::InvalidConfiguration(format!(
                "Venturi {name} must lie in [1, 45] degrees"
            )))
        }
    };

    let inlet_width = validate_positive("inlet width", metadata.inlet_width_m)?;
    let throat_width = validate_positive("throat width", metadata.throat_width_m)?;
    let throat_height = validate_positive("throat height", metadata.throat_height_m)?;
    let outlet_width = validate_positive("outlet width", metadata.outlet_width_m)?;
    let throat_length = validate_positive("throat length", metadata.throat_length_m)?;
    let convergent_angle =
        validate_angle("convergent half-angle", metadata.convergent_half_angle_deg)?;
    let divergent_angle =
        validate_angle("divergent half-angle", metadata.divergent_half_angle_deg)?;

    let hydraulic_diameter =
        |width_m: f64, height_m: f64| 2.0 * width_m * height_m / (width_m + height_m);
    let inlet_d = hydraulic_diameter(inlet_width, throat_height);
    let throat_d = hydraulic_diameter(throat_width, throat_height);
    let outlet_d = hydraulic_diameter(outlet_width, throat_height);
    let conv_len = venturi_taper_length_m(inlet_width, throat_width, convergent_angle)?;
    let diff_len = venturi_taper_length_m(outlet_width, throat_width, divergent_angle)?;
    let total_length = throat_length + conv_len + diff_len;

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
    let throat_len_t = T::from_f64(throat_length).ok_or_else(|| {
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
                convergent_angle,
            )?,
        })
        .with_expansion(ExpansionType::Gradual {
            half_angle_deg: divergent_angle,
        });
    model.throat_roughness =
        T::from_f64(1e-7).expect("Mathematical constant conversion compromised");

    let conds = ModelFlowConditions::from_flow_rate(
        T::from_f64(1e-12).expect("Mathematical constant conversion compromised"),
    );
    model.calculate_coefficients(fluid, &conds)
}

#[cfg(test)]
mod tests {
    use super::venturi_coefficients;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::database::water_20c;
    use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;

    #[test]
    fn venturi_coefficients_matches_manual_model() -> cfd_core::error::Result<()> {
        use crate::physics::resistance::models::{
            ExpansionType, FlowConditions as ModelFlowConditions, ResistanceModel, VenturiGeometry,
            VenturiModel,
        };
        use crate::{discharge_coefficient_from_convergent_half_angle_deg, venturi_taper_length_m};

        let metadata = VenturiGeometryMetadata {
            throat_width_m: 4.0e-4,
            throat_height_m: 2.0e-4,
            throat_length_m: 8.0e-4,
            inlet_width_m: 1.0e-3,
            outlet_width_m: 1.2e-3,
            convergent_half_angle_deg: 7.0,
            divergent_half_angle_deg: 5.0,
            throat_position: 0.5,
        };
        let fluid = water_20c::<f64>()?;

        let (r, k) = venturi_coefficients(&metadata, &fluid)?;

        let hydraulic_diameter =
            |width_m: f64, height_m: f64| 2.0 * width_m * height_m / (width_m + height_m);
        let inlet_d = hydraulic_diameter(metadata.inlet_width_m, metadata.throat_height_m);
        let throat_d = hydraulic_diameter(metadata.throat_width_m, metadata.throat_height_m);
        let outlet_d = hydraulic_diameter(metadata.outlet_width_m, metadata.throat_height_m);
        let conv_len = venturi_taper_length_m(
            metadata.inlet_width_m,
            metadata.throat_width_m,
            metadata.convergent_half_angle_deg,
        )?;
        let diff_len = venturi_taper_length_m(
            metadata.outlet_width_m,
            metadata.throat_width_m,
            metadata.divergent_half_angle_deg,
        )?;
        let total_length = metadata.throat_length_m + conv_len + diff_len;

        let inlet_d_t = f64::from(inlet_d);
        let throat_d_t = f64::from(throat_d);
        let outlet_d_t = f64::from(outlet_d);
        let throat_len_t = f64::from(metadata.throat_length_m);
        let total_len_t = f64::from(total_length);

        let manual =
            VenturiModel::new(inlet_d_t, throat_d_t, outlet_d_t, throat_len_t, total_len_t)
                .with_geometry(VenturiGeometry::Custom {
                    discharge_coefficient: discharge_coefficient_from_convergent_half_angle_deg(
                        metadata.convergent_half_angle_deg,
                    )?,
                })
                .with_expansion(ExpansionType::Gradual {
                    half_angle_deg: metadata.divergent_half_angle_deg,
                });

        let conditions = ModelFlowConditions::from_flow_rate(1.0e-12);
        let (manual_r, manual_k) = manual.calculate_coefficients(&fluid, &conditions)?;

        assert_relative_eq!(r, manual_r, epsilon = manual_r.abs().max(1.0) * 1e-12);
        assert_relative_eq!(k, manual_k, epsilon = manual_k.abs().max(1.0) * 1e-12);

        Ok(())
    }

    #[test]
    fn venturi_coefficients_avoid_zero_velocity_fallback() -> cfd_core::error::Result<()> {
        use crate::physics::resistance::models::{
            ExpansionType, FlowConditions as ModelFlowConditions, ResistanceModel, VenturiGeometry,
            VenturiModel,
        };
        use crate::{discharge_coefficient_from_convergent_half_angle_deg, venturi_taper_length_m};

        let metadata = VenturiGeometryMetadata {
            throat_width_m: 4.0e-4,
            throat_height_m: 2.0e-4,
            throat_length_m: 8.0e-4,
            inlet_width_m: 1.0e-3,
            outlet_width_m: 1.2e-3,
            convergent_half_angle_deg: 12.0,
            divergent_half_angle_deg: 7.0,
            throat_position: 0.5,
        };
        let fluid = water_20c::<f64>()?;

        let (_, builder_k) = venturi_coefficients(&metadata, &fluid)?;

        let hydraulic_diameter =
            |width_m: f64, height_m: f64| 2.0 * width_m * height_m / (width_m + height_m);
        let inlet_d = hydraulic_diameter(metadata.inlet_width_m, metadata.throat_height_m);
        let throat_d = hydraulic_diameter(metadata.throat_width_m, metadata.throat_height_m);
        let outlet_d = hydraulic_diameter(metadata.outlet_width_m, metadata.throat_height_m);
        let conv_len = venturi_taper_length_m(
            metadata.inlet_width_m,
            metadata.throat_width_m,
            metadata.convergent_half_angle_deg,
        )?;
        let diff_len = venturi_taper_length_m(
            metadata.outlet_width_m,
            metadata.throat_width_m,
            metadata.divergent_half_angle_deg,
        )?;
        let total_length = metadata.throat_length_m + conv_len + diff_len;

        let manual = VenturiModel::new(
            inlet_d,
            throat_d,
            outlet_d,
            metadata.throat_length_m,
            total_length,
        )
        .with_geometry(VenturiGeometry::Custom {
            discharge_coefficient: discharge_coefficient_from_convergent_half_angle_deg(
                metadata.convergent_half_angle_deg,
            )?,
        })
        .with_expansion(ExpansionType::Gradual {
            half_angle_deg: metadata.divergent_half_angle_deg,
        });

        let flow_probe = ModelFlowConditions::from_flow_rate(1.0e-12);
        let (_, flow_probe_k) = manual.calculate_coefficients(&fluid, &flow_probe)?;

        let mut zero_velocity = ModelFlowConditions::new(0.0);
        zero_velocity.flow_rate = Some(1.0e-12);
        let (_, zero_velocity_k) = manual.calculate_coefficients(&fluid, &zero_velocity)?;

        assert_relative_eq!(
            builder_k,
            flow_probe_k,
            epsilon = flow_probe_k.abs().max(1.0) * 1e-12
        );
        assert!(builder_k > 0.0);
        assert_eq!(zero_velocity_k, 0.0);
        assert!(builder_k > zero_velocity_k);

        Ok(())
    }
}
