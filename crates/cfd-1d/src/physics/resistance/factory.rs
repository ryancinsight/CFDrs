//! Factory for creating resistance models.

use super::models::{
    BendType, DarcyWeisbachModel, ExpansionType, HagenPoiseuilleModel, RectangularChannelModel,
    SerpentineCrossSection, SerpentineModel, VenturiGeometry, VenturiModel,
};
use super::traits::{scalar_to_f64, ResistanceScalar};

/// Resistance model factory for creating standard models
pub struct ResistanceModelFactory;

impl ResistanceModelFactory {
    /// Create Hagen-Poiseuille model for circular channel
    pub fn hagen_poiseuille<T: ResistanceScalar>(
        diameter: T,
        length: T,
    ) -> HagenPoiseuilleModel<T> {
        HagenPoiseuilleModel { diameter, length }
    }

    /// Create rectangular channel model
    pub fn rectangular_channel<T: ResistanceScalar>(
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
    pub fn darcy_weisbach<T: ResistanceScalar>(
        hydraulic_diameter: T,
        length: T,
        roughness: T,
    ) -> DarcyWeisbachModel<T> {
        DarcyWeisbachModel::circular(hydraulic_diameter, length, roughness)
    }

    /// Create Darcy-Weisbach model for turbulent flow in a circular channel
    pub fn darcy_weisbach_circular<T: ResistanceScalar>(
        diameter: T,
        length: T,
        roughness: T,
    ) -> DarcyWeisbachModel<T> {
        DarcyWeisbachModel::circular(diameter, length, roughness)
    }

    /// Create serpentine channel resistance model with circular cross-section.
    ///
    /// Accounts for Dean flow corrections in curved channels and bend losses.
    /// Uses White (1929) / Ito (1959) curvature enhancement correlations and
    /// Idelchik (2007) bend loss coefficients.
    ///
    /// # Arguments
    /// - `diameter`: Channel diameter \[m]
    /// - `straight_length`: Total length of all straight segments \[m]
    /// - `num_segments`: Number of straight segments (bends = segments - 1)
    /// - `bend_radius`: Radius of curvature of bends \[m]
    pub fn serpentine_circular<T: ResistanceScalar>(
        diameter: T,
        straight_length: T,
        num_segments: usize,
        bend_radius: T,
    ) -> SerpentineModel<T> {
        let dh_f64 = scalar_to_f64::<T>(diameter);
        let br_f64 = scalar_to_f64::<T>(bend_radius);
        let ratio = if dh_f64 > 0.0 { br_f64 / dh_f64 } else { 2.0 };

        SerpentineModel {
            straight_length,
            num_segments,
            cross_section: SerpentineCrossSection::Circular { diameter: dh_f64 },
            bend_radius,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: ratio,
            },
        }
    }

    /// Create serpentine channel resistance model with rectangular cross-section.
    ///
    /// Applies Shah-London f·Re corrections for rectangular ducts and
    /// Dean flow curvature enhancement.
    ///
    /// # Arguments
    /// - `width`: Channel width \[m]
    /// - `height`: Channel height (depth) \[m]
    /// - `straight_length`: Total length of all straight segments \[m]
    /// - `num_segments`: Number of straight segments (bends = segments - 1)
    /// - `bend_radius`: Radius of curvature of bends \[m]
    pub fn serpentine_rectangular<T: ResistanceScalar>(
        width: T,
        height: T,
        straight_length: T,
        num_segments: usize,
        bend_radius: T,
    ) -> SerpentineModel<T> {
        let w = scalar_to_f64::<T>(width);
        let h = scalar_to_f64::<T>(height);
        let dh = 2.0 * w * h / (w + h);
        let br = scalar_to_f64::<T>(bend_radius);
        let ratio = if dh > 0.0 { br / dh } else { 2.0 };

        SerpentineModel {
            straight_length,
            num_segments,
            cross_section: SerpentineCrossSection::Rectangular {
                width: w,
                height: h,
            },
            bend_radius,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: ratio,
            },
        }
    }

    /// Create symmetric Venturi tube model (outlet diameter = inlet diameter).
    ///
    /// Uses ISO 5167-4 discharge coefficient correlations and Borda-Carnot
    /// expansion loss model.
    ///
    /// # Arguments
    /// - `inlet_diameter`: Upstream (inlet) pipe diameter \[m]
    /// - `throat_diameter`: Throat (constriction) diameter \[m]
    /// - `throat_length`: Length of the throat section \[m]
    /// - `total_length`: Total device length \[m]
    pub fn venturi_symmetric<T: ResistanceScalar>(
        inlet_diameter: T,
        throat_diameter: T,
        throat_length: T,
        total_length: T,
    ) -> VenturiModel<T> {
        VenturiModel::symmetric(inlet_diameter, throat_diameter, throat_length, total_length)
    }

    /// Create millifluidic Venturi with typical parameters for blood flow.
    ///
    /// Pre-configured with:
    /// - Custom discharge coefficient (0.97, typical for machined millifluidic devices)
    /// - Gradual 5° half-angle diffuser
    /// - Smooth throat (zero roughness)
    ///
    /// # Arguments
    /// - `inlet_diameter`: Upstream pipe diameter \[m]
    /// - `throat_diameter`: Throat diameter \[m]
    /// - `throat_length`: Length of the throat section \[m]
    pub fn venturi_millifluidic<T: ResistanceScalar>(
        inlet_diameter: T,
        throat_diameter: T,
        throat_length: T,
    ) -> VenturiModel<T> {
        VenturiModel::millifluidic(inlet_diameter, throat_diameter, throat_length)
    }

    /// Create Venturi model with full configuration.
    ///
    /// # Arguments
    /// - `inlet_diameter`: Upstream pipe diameter \[m]
    /// - `throat_diameter`: Throat diameter \[m]
    /// - `outlet_diameter`: Downstream pipe diameter \[m]
    /// - `throat_length`: Length of the throat section \[m]
    /// - `total_length`: Total device length \[m]
    /// - `geometry`: Venturi geometry type (determines discharge coefficient)
    /// - `expansion`: Expansion type (determines recovery efficiency)
    pub fn venturi_custom<T: ResistanceScalar>(
        inlet_diameter: T,
        throat_diameter: T,
        outlet_diameter: T,
        throat_length: T,
        total_length: T,
        geometry: VenturiGeometry,
        expansion: ExpansionType,
    ) -> VenturiModel<T> {
        VenturiModel::new(
            inlet_diameter,
            throat_diameter,
            outlet_diameter,
            throat_length,
            total_length,
        )
        .with_geometry(geometry)
        .with_expansion(expansion)
    }
}
