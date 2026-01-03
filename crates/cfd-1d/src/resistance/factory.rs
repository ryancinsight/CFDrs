//! Factory for creating resistance models.

use super::models::{DarcyWeisbachModel, HagenPoiseuilleModel, RectangularChannelModel};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Resistance model factory for creating standard models
pub struct ResistanceModelFactory;

impl ResistanceModelFactory {
    /// Create Hagen-Poiseuille model for circular channel
    pub fn hagen_poiseuille<T: RealField + Copy + FromPrimitive + Copy>(
        diameter: T,
        length: T,
    ) -> HagenPoiseuilleModel<T> {
        HagenPoiseuilleModel { diameter, length }
    }

    /// Create rectangular channel model
    pub fn rectangular_channel<T: RealField + Copy + FromPrimitive + Copy>(
        width: T,
        height: T,
        length: T,
    ) -> RectangularChannelModel<T> {
        RectangularChannelModel {
            width,
            height,
            length,
        }
    }

    /// Create Darcy-Weisbach model for turbulent flow in any geometry
    pub fn darcy_weisbach<T: RealField + Copy + FromPrimitive>(
        hydraulic_diameter: T,
        length: T,
        roughness: T,
    ) -> DarcyWeisbachModel<T> {
           DarcyWeisbachModel::circular(hydraulic_diameter, length, roughness)
    }

    /// Create Darcy-Weisbach model for turbulent flow in a circular channel
    pub fn darcy_weisbach_circular<T: RealField + Copy + FromPrimitive>(
        diameter: T,
        length: T,
        roughness: T,
    ) -> DarcyWeisbachModel<T> {
        DarcyWeisbachModel::circular(diameter, length, roughness)
    }
}
