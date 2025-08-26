//! Resistance calculator with model selection and validation

use super::models::{DarcyWeisbachModel, HagenPoiseuilleModel, RectangularChannelModel};
use super::{ChannelGeometry, FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
/// Resistance calculator with model selection and validation
pub struct ResistanceCalculator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}
impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceCalculator<T> {
    /// Create new resistance calculator
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
    /// Calculate resistance with automatic model selection
    pub fn calculate_auto(
        &self,
        geometry: &ChannelGeometry<T>,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        match geometry {
            ChannelGeometry::Circular { diameter, length } => {
                let model = HagenPoiseuilleModel {
                    diameter: *diameter,
                    length: *length,
                };
                model.calculate_resistance(fluid, conditions)
            }
            ChannelGeometry::Rectangular {
                width,
                height,
                length,
            } => {
                let model = RectangularChannelModel {
                    width: *width,
                    height: *height,
    /// Calculate resistance using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille(
        diameter: T,
        length: T,
        let model = HagenPoiseuilleModel { diameter, length };
        model.calculate_resistance(fluid, conditions)
    /// Calculate resistance using rectangular channel model
    pub fn calculate_rectangular(
        width: T,
        height: T,
        let model = RectangularChannelModel {
            width,
            height,
            length,
        };
    /// Calculate resistance using Darcy-Weisbach model
    pub fn calculate_darcy_weisbach(
        hydraulic_diameter: T,
        roughness: T,
        let model = DarcyWeisbachModel {
            hydraulic_diameter,
            roughness,
impl<T: RealField + Copy + FromPrimitive + num_traits::Float> Default for ResistanceCalculator<T> {
    fn default() -> Self {
        Self::new()
