//! Factory for creating resistance models

use super::models::{DarcyWeisbachModel, HagenPoiseuilleModel, RectangularChannelModel};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
/// Methods for combining multiple resistance models
#[derive(Debug, Clone, PartialEq)]
pub enum CombinationMethod {
    /// Add resistances in series
    Series,
    /// Add resistances in parallel
    Parallel,
    /// Custom weighted combination
    Weighted(Vec<f64>),
}
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
    ) -> RectangularChannelModel<T> {
        RectangularChannelModel {
            width,
            height,
            length,
        }
    /// Create Darcy-Weisbach model for turbulent flow
    pub fn darcy_weisbach<T: RealField + Copy + FromPrimitive + Copy>(
        hydraulic_diameter: T,
        roughness: T,
    ) -> DarcyWeisbachModel<T> {
        DarcyWeisbachModel {
            hydraulic_diameter,
            roughness,

    }


}
}
}
