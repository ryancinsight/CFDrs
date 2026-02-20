use crate::resistance::geometry::ChannelGeometry;
use crate::resistance::models::{
    DarcyWeisbachModel, FlowConditions, HagenPoiseuilleModel, MembranePoreModel,
    RectangularChannelModel, ResistanceModel, SerpentineModel, VenturiModel,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Calculate resistance with automatic model selection
pub fn calculate_auto<T, F>(
    geometry: &ChannelGeometry<T>,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let mut local_conditions = conditions.clone();

    // Compute Reynolds number if not provided.
    if local_conditions.reynolds_number.is_none() {
        let state = fluid.properties_at(local_conditions.temperature, local_conditions.pressure)?;
        let density = state.density;
        let viscosity = state.dynamic_viscosity;

        let velocity = if let Some(v) = local_conditions.velocity {
            v
        } else if let Some(q) = local_conditions.flow_rate {
            let area = geometry.cross_sectional_area()?;
            if area <= T::zero() {
                return Err(Error::InvalidConfiguration(
                    "Channel area must be positive to compute Reynolds number".to_string(),
                ));
            }
            q / area
        } else {
            return Err(Error::InvalidConfiguration(
                "Automatic resistance selection requires either velocity or flow_rate (or an explicit Reynolds number)".to_string(),
            ));
        };

        let dh = geometry.hydraulic_diameter()?;
        let re = density * velocity * dh / viscosity;

        local_conditions.velocity = Some(velocity);
        local_conditions.reynolds_number = Some(re);
    }

    match geometry {
        ChannelGeometry::Circular { diameter, length } => {
            let hp = HagenPoiseuilleModel {
                diameter: *diameter,
                length: *length,
            };
            if hp.is_applicable(&local_conditions) {
                hp.validate_invariants(fluid, &local_conditions)?;
                return hp.calculate_resistance(fluid, &local_conditions);
            }

            // Default to a smooth pipe if roughness is not specified by the geometry.
            let dw = DarcyWeisbachModel::circular(*diameter, *length, T::zero());
            if dw.is_applicable(&local_conditions) {
                dw.validate_invariants(fluid, &local_conditions)?;
                return dw.calculate_resistance(fluid, &local_conditions);
            }

            Err(Error::InvalidConfiguration(
                "No applicable resistance model for circular channel at the given Reynolds number".to_string(),
            ))
        }
        ChannelGeometry::Rectangular {
            width,
            height,
            length,
        } => {
            let rect = RectangularChannelModel {
                width: *width,
                height: *height,
                length: *length,
            };
            if rect.is_applicable(&local_conditions) {
                rect.validate_invariants(fluid, &local_conditions)?;
                rect.calculate_resistance(fluid, &local_conditions)
            } else {
                Err(Error::InvalidConfiguration(
                    "Rectangular-channel resistance currently supports laminar flow only; provide a laminar Reynolds number or use a different model".to_string(),
                ))
            }
        }
        _ => Err(Error::InvalidConfiguration(
            "No resistance model available for this geometry type; use Circular or Rectangular".to_string(),
        )),
    }
}

pub fn calculate_hagen_poiseuille<T, F>(
    diameter: T,
    length: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = HagenPoiseuilleModel { diameter, length };
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_rectangular<T, F>(
    width: T,
    height: T,
    length: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = RectangularChannelModel {
        width,
        height,
        length,
    };
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_darcy_weisbach<T, F>(
    hydraulic_diameter: T,
    length: T,
    roughness: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = DarcyWeisbachModel::circular(hydraulic_diameter, length, roughness);
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_serpentine_circular<T, F>(
    diameter: T,
    straight_length: T,
    num_segments: usize,
    bend_radius: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = crate::resistance::factory::ResistanceModelFactory::serpentine_circular(
        diameter,
        straight_length,
        num_segments,
        bend_radius,
    );
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_serpentine_rectangular<T, F>(
    width: T,
    height: T,
    straight_length: T,
    num_segments: usize,
    bend_radius: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = crate::resistance::factory::ResistanceModelFactory::serpentine_rectangular(
        width,
        height,
        straight_length,
        num_segments,
        bend_radius,
    );
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_venturi<T, F>(
    inlet_diameter: T,
    throat_diameter: T,
    throat_length: T,
    total_length: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = VenturiModel::symmetric(
        inlet_diameter,
        throat_diameter,
        throat_length,
        total_length,
    );
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

pub fn calculate_membrane_porous<T, F>(
    thickness: T,
    width: T,
    height: T,
    pore_radius: T,
    porosity: T,
    fluid: &F,
    conditions: &FlowConditions<T>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    let model = MembranePoreModel::new(thickness, width, height, pore_radius, porosity);
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}
