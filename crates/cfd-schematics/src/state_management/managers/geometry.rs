//! Geometry parameter manager.

use crate::config::constants::primitives as constants;
use crate::state_management::{
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
    parameters::{ConfigurableParameter, ParameterMetadata},
};

use super::ParameterManager;

/// Parameter manager for geometry parameters
#[derive(Debug)]
pub struct GeometryParameterManager {
    /// Wall clearance parameter
    wall_clearance: ConfigurableParameter<f64>,

    /// Channel width parameter
    channel_width: ConfigurableParameter<f64>,

    /// Channel height parameter
    channel_height: ConfigurableParameter<f64>,
}

impl GeometryParameterManager {
    /// Create a new geometry parameter manager
    #[must_use]
    pub fn new() -> Self {
        let wall_clearance = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(
                    constants::MIN_WALL_CLEARANCE,
                    constants::MAX_WALL_CLEARANCE,
                ),
            ]),
            ParameterMetadata::new(
                "wall_clearance",
                "Minimum distance between channels and walls",
                "geometry_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let channel_width = ConfigurableParameter::new(
            1.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(
                    constants::MIN_CHANNEL_WIDTH,
                    constants::MAX_CHANNEL_WIDTH,
                ),
            ]),
            ParameterMetadata::new(
                "channel_width",
                "Width of microfluidic channels",
                "geometry_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let channel_height = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(
                    constants::MIN_CHANNEL_HEIGHT,
                    constants::MAX_CHANNEL_HEIGHT,
                ),
            ]),
            ParameterMetadata::new(
                "channel_height",
                "Height of microfluidic channels",
                "geometry_parameters",
            )
            .with_units("mm"),
        );

        Self {
            wall_clearance,
            channel_width,
            channel_height,
        }
    }
}

impl Default for GeometryParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for GeometryParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "wall_clearance" => Ok(Box::new(*self.wall_clearance.get_raw_value())),
            "channel_width" => Ok(Box::new(*self.channel_width.get_raw_value())),
            "channel_height" => Ok(Box::new(*self.channel_height.get_raw_value())),
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "wall_clearance" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.wall_clearance.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "channel_width" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.channel_width.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "channel_height" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.channel_height.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "wall_clearance".to_string(),
            "channel_width".to_string(),
            "channel_height".to_string(),
        ]
    }

    fn has_parameter(&self, name: &str) -> bool {
        self.parameter_names().contains(&name.to_string())
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.wall_clearance.validate()?;
        self.channel_width.validate()?;
        self.channel_height.validate()?;

        let wall_clearance = *self.wall_clearance.get_raw_value();
        let channel_width = *self.channel_width.get_raw_value();
        if wall_clearance >= channel_width {
            return Err(ParameterError::invalid_value(
                "wall_clearance",
                &wall_clearance,
                "must be less than channel_width",
            ));
        }

        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "wall_clearance" => Ok(self.wall_clearance.metadata()),
            "channel_width" => Ok(self.channel_width.metadata()),
            "channel_height" => Ok(self.channel_height.metadata()),
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "geometry"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.wall_clearance.reset(reason)?;
        self.channel_width.reset(reason)?;
        self.channel_height.reset(reason)?;
        Ok(())
    }
}
