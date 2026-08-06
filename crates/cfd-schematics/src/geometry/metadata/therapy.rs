use aequitas::systems::si::quantities::{
    Angle, Dimensionless, Length, MassDensity, Pressure, Velocity, VolumetricFlowRate,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiGeometryMetadata {
    pub throat_width_m: Length<f64>,
    pub throat_height_m: Length<f64>,
    pub throat_length_m: Length<f64>,
    pub inlet_width_m: Length<f64>,
    pub outlet_width_m: Length<f64>,
    pub convergent_half_angle: Angle<f64>,
    pub divergent_half_angle: Angle<f64>,
    #[serde(default = "default_throat_position")]
    pub throat_position: Dimensionless<f64>,
}

fn default_throat_position() -> Dimensionless<f64> {
    Dimensionless::from_base(0.5)
}

crate::impl_metadata!(VenturiGeometryMetadata, "VenturiGeometryMetadata");

#[derive(Debug, Clone, PartialEq)]
pub struct CascadeParams {
    pub n_levels: u8,
    pub center_frac: f64,
}

crate::impl_metadata!(CascadeParams, "CascadeParams");

#[derive(Debug, Clone, PartialEq)]
pub struct IncrementalFiltrationParams {
    pub n_pretri: u8,
    pub pretri_center_frac: f64,
    pub terminal_tri_center_frac: f64,
    pub bi_treat_frac: f64,
    pub outlet_tail_length_m: Length<f64>,
}

crate::impl_metadata!(IncrementalFiltrationParams, "IncrementalFiltrationParams");

#[derive(Debug, Clone, PartialEq)]
pub struct AsymmetricTrifurcationParams {
    pub center_frac: f64,
    pub left_frac: f64,
    pub right_frac: f64,
}

crate::impl_metadata!(AsymmetricTrifurcationParams, "AsymmetricTrifurcationParams");

#[derive(Debug, Clone, PartialEq)]
pub struct ChannelVenturiSpec {
    pub n_throats: u8,
    pub is_ctc_stream: bool,
    pub throat_width_m: Length<f64>,
    pub height_m: Length<f64>,
    pub inter_throat_spacing_m: Length<f64>,
}

impl ChannelVenturiSpec {
    #[must_use]
    pub fn cumulative_cavitation_dose(
        &self,
        cav_potential: Dimensionless<f64>,
    ) -> Dimensionless<f64> {
        let p = cav_potential.into_base().clamp(0.0, 1.0);
        if p <= 0.0 || self.n_throats == 0 {
            return Dimensionless::from_base(0.0);
        }
        Dimensionless::from_base((1.0 - (1.0 - p).powi(i32::from(self.n_throats))).clamp(0.0, 1.0))
    }

    #[must_use]
    pub fn total_throat_pressure_drop(
        &self,
        flow_rate: VolumetricFlowRate<f64>,
        blood_density: MassDensity<f64>,
        diffuser_coeff: Dimensionless<f64>,
    ) -> Pressure<f64> {
        if self.n_throats == 0 {
            return Pressure::from_base(0.0);
        }
        let throat_width = self.throat_width_m.into_base();
        let height = self.height_m.into_base();
        let flow_rate = flow_rate.into_base();
        let blood_density = blood_density.into_base();
        let diffuser_coeff = diffuser_coeff.into_base();
        let inlet_area = (throat_width * 2.0) * height;
        let throat_area = throat_width * height;
        let v_in = flow_rate / inlet_area.max(1e-18);
        let v_throat = flow_rate / throat_area.max(1e-18);
        let single_dp = 0.5
            * blood_density
            * (v_throat * v_throat - v_in * v_in).max(0.0)
            * (1.0 - diffuser_coeff);
        Pressure::from_base(single_dp * f64::from(self.n_throats))
    }
}

crate::impl_metadata!(ChannelVenturiSpec, "ChannelVenturiSpec");

#[derive(Debug, Clone, PartialEq)]
pub struct FdaCavitationCompliance {
    pub mi_equiv: Dimensionless<f64>,
    pub fda_mi_compliant: bool,
    pub therapeutic_mi_lower: Dimensionless<f64>,
    pub in_therapeutic_window: bool,
}

impl FdaCavitationCompliance {
    pub const FDA_MI_LIMIT: Dimensionless<f64> = Dimensionless::from_base(1.9);
    pub const INERTIAL_CAV_THRESHOLD_MI: Dimensionless<f64> = Dimensionless::from_base(0.3);
    pub const SOUND_SPEED_BLOOD_M_S: Velocity<f64> = Velocity::from_base(1540.0);

    #[must_use]
    pub fn from_throat_dp(dp_throat: Pressure<f64>, blood_density: MassDensity<f64>) -> Self {
        let dp_throat = dp_throat.into_base();
        let blood_density = blood_density.into_base();
        let v_rms = (2.0 * dp_throat.max(0.0) / blood_density.max(1.0)).sqrt();
        let mi_equiv = Dimensionless::from_base(v_rms / Self::SOUND_SPEED_BLOOD_M_S.into_base());
        Self {
            mi_equiv,
            fda_mi_compliant: mi_equiv.into_base() < Self::FDA_MI_LIMIT.into_base(),
            therapeutic_mi_lower: Self::INERTIAL_CAV_THRESHOLD_MI,
            in_therapeutic_window: (Self::INERTIAL_CAV_THRESHOLD_MI.into_base()
                ..Self::FDA_MI_LIMIT.into_base())
                .contains(&mi_equiv.into_base()),
        }
    }
}

crate::impl_metadata!(FdaCavitationCompliance, "FdaCavitationCompliance");

#[cfg(test)]
mod tests {
    use super::ChannelVenturiSpec;
    use aequitas::systems::si::quantities::{
        Dimensionless, Length, MassDensity, Pressure, VolumetricFlowRate,
    };

    #[test]
    fn cumulative_cavitation_dose_compounds_per_throat() {
        let spec = ChannelVenturiSpec {
            n_throats: 2,
            is_ctc_stream: true,
            throat_width_m: Length::from_base(0.5),
            height_m: Length::from_base(1.0),
            inter_throat_spacing_m: Length::from_base(1.0),
        };

        assert_eq!(
            spec.cumulative_cavitation_dose(Dimensionless::from_base(0.5))
                .into_base(),
            0.75
        );
    }

    #[test]
    fn total_throat_pressure_drop_preserves_pressure_dimension() {
        let spec = ChannelVenturiSpec {
            n_throats: 2,
            is_ctc_stream: false,
            throat_width_m: Length::from_base(0.5),
            height_m: Length::from_base(1.0),
            inter_throat_spacing_m: Length::from_base(1.0),
        };

        let pressure = spec.total_throat_pressure_drop(
            VolumetricFlowRate::from_base(1.0),
            MassDensity::from_base(2.0),
            Dimensionless::from_base(0.25),
        );

        assert_eq!(pressure.into_base(), 4.5);
    }

    #[test]
    fn fda_cavitation_compliance_accepts_typed_pressure_and_density() {
        let compliance = super::FdaCavitationCompliance::from_throat_dp(
            Pressure::from_base(2.0),
            MassDensity::from_base(2.0),
        );

        assert_eq!(compliance.mi_equiv.into_base(), 2.0_f64.sqrt() / 1540.0);
        assert_eq!(
            compliance.therapeutic_mi_lower,
            Dimensionless::from_base(0.3)
        );
        assert!(compliance.fda_mi_compliant);
        assert!(!compliance.in_therapeutic_window);
    }
}
