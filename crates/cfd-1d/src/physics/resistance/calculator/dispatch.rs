//! Resistance model dispatch functions.
//!
//! Provides convenience free functions that construct and call the appropriate
//! `ResistanceModel` implementation for each geometry, as well as an
//! automatic model-selection function (`calculate_auto`) that chooses between
//! validated geometry-owned models and requires explicit opt-in when a caller
//! wants a different surrogate such as Darcy-Weisbach.
use crate::physics::resistance::geometry::ChannelGeometry;
use crate::physics::resistance::models::{
    DarcyWeisbachModel, FlowConditions, HagenPoiseuilleModel, MembranePoreModel,
    RectangularChannelModel, ResistanceModel, VenturiModel,
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
    super::populate_shear_aware_conditions(geometry, fluid, &mut local_conditions)?;

    match geometry {
        ChannelGeometry::Circular { diameter, length } => {
            let hp = HagenPoiseuilleModel {
                diameter: *diameter,
                length: *length,
            };
            if hp.is_applicable(&local_conditions) {
                return hp.calculate_resistance(fluid, &local_conditions);
            }

            // Default to a smooth pipe if roughness is not specified by the geometry.
            let dw = DarcyWeisbachModel::circular(*diameter, *length, T::zero());
            if dw.is_applicable(&local_conditions) {
                return dw.calculate_resistance(fluid, &local_conditions);
            }

            Err(Error::InvalidConfiguration(
                "No applicable resistance model for circular channel at the given Reynolds number"
                    .to_string(),
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
                return rect.calculate_resistance(fluid, &local_conditions);
            }

            Err(super::rectangular_auto_selection_error(&local_conditions))
        }
        _ => Err(Error::InvalidConfiguration(
            "No resistance model available for this geometry type; use Circular or Rectangular"
                .to_string(),
        )),
    }
}

/// Calculate hydraulic resistance using the Hagen-Poiseuille model for a circular pipe.
///
/// # Errors
/// Returns `Err` if `diameter ≤ 0`, `length ≤ 0`, or if the fluid properties are invalid.
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

/// Calculate hydraulic resistance for a rectangular cross-section channel (laminar flow).
///
/// Uses the infinite-series Shah & London correction for the aspect ratio.
/// Valid only for laminar flow (Re < 2300).
///
/// # Errors
/// Returns `Err` if any dimension is non-positive or if Re ≥ 2300.
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

/// Calculate hydraulic resistance using the Darcy-Weisbach model (circular section).
///
/// Uses the Colebrook-White / Haaland equation for the friction factor.
/// Applies to both laminar and turbulent flow regimes.
///
/// # Errors
/// Returns `Err` if `hydraulic_diameter ≤ 0`, `length ≤ 0`, or `roughness < 0`.
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

/// Calculate resistance for a circular serpentine channel including Dean flow correction.
///
/// Accounts for secondary Dean vortices in curved segments using the Ito (1959) correlation.
///
/// # Errors
/// Returns `Err` if any input dimension is non-positive or `num_segments == 0`.
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
    let model = crate::physics::resistance::factory::ResistanceModelFactory::serpentine_circular(
        diameter,
        straight_length,
        num_segments,
        bend_radius,
    );
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

/// Calculate resistance for a rectangular serpentine channel including Dean flow correction.
///
/// Uses the rectangular hydraulic diameter `D_h = 2wh/(w+h)` for the Dean number.
///
/// # Errors
/// Returns `Err` if any dimension is non-positive or `num_segments == 0`.
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
    let model = crate::physics::resistance::factory::ResistanceModelFactory::serpentine_rectangular(
        width,
        height,
        straight_length,
        num_segments,
        bend_radius,
    );
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

/// Calculate total hydraulic resistance for a Venturi (converging-diverging) tube.
///
/// Includes: convergent pressure recovery, throat viscous friction, and Borda-Carnot
/// expansion loss in the divergent section (ISO 5167-4 discharge coefficients).
///
/// # Errors
/// Returns `Err` if `throat_diameter ≥ inlet_diameter`, or any length is non-positive.
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
    let model =
        VenturiModel::symmetric(inlet_diameter, throat_diameter, throat_length, total_length);
    model.validate_invariants(fluid, conditions)?;
    model.calculate_resistance(fluid, conditions)
}

/// Calculate hydraulic resistance through a porous membrane modelled as parallel cylindrical pores.
///
/// Uses the Hagen-Poiseuille law per pore: `R_total = R_pore / N_pores`, where
/// `N_pores = porosity × A_membrane / (π r_pore²)`.
///
/// # Errors
/// Returns `Err` if `pore_radius ≤ 0`, `porosity ∉ (0, 1]`, or any dimension is non-positive.
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
