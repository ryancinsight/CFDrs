use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiGeometryMetadata {
    pub throat_width_m: f64,
    pub throat_height_m: f64,
    pub throat_length_m: f64,
    pub inlet_width_m: f64,
    pub outlet_width_m: f64,
    pub convergent_half_angle_deg: f64,
    pub divergent_half_angle_deg: f64,
    #[serde(default = "default_throat_position")]
    pub throat_position: f64,
}

fn default_throat_position() -> f64 {
    0.5
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
    pub outlet_tail_length_m: f64,
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
    pub throat_width_m: f64,
    pub height_m: f64,
    pub inter_throat_spacing_m: f64,
}

impl ChannelVenturiSpec {
    #[must_use]
    pub fn cumulative_cavitation_dose(&self, cav_potential: f64) -> f64 {
        let p = cav_potential.clamp(0.0, 1.0);
        if p <= 0.0 || self.n_throats == 0 {
            return 0.0;
        }
        (1.0 - (1.0 - p).powi(i32::from(self.n_throats))).clamp(0.0, 1.0)
    }

    #[must_use]
    pub fn total_throat_pressure_drop_pa(
        &self,
        flow_m3_s: f64,
        blood_density_kg_m3: f64,
        diffuser_coeff: f64,
    ) -> f64 {
        if self.n_throats == 0 {
            return 0.0;
        }
        let inlet_area = (self.throat_width_m * 2.0) * self.height_m;
        let throat_area = self.throat_width_m * self.height_m;
        let v_in = flow_m3_s / inlet_area.max(1e-18);
        let v_throat = flow_m3_s / throat_area.max(1e-18);
        let single_dp = 0.5
            * blood_density_kg_m3
            * (v_throat * v_throat - v_in * v_in).max(0.0)
            * (1.0 - diffuser_coeff);
        single_dp * f64::from(self.n_throats)
    }
}

crate::impl_metadata!(ChannelVenturiSpec, "ChannelVenturiSpec");

#[derive(Debug, Clone, PartialEq)]
pub struct FdaCavitationCompliance {
    pub mi_equiv: f64,
    pub fda_mi_compliant: bool,
    pub therapeutic_mi_lower: f64,
    pub in_therapeutic_window: bool,
}

impl FdaCavitationCompliance {
    pub const FDA_MI_LIMIT: f64 = 1.9;
    pub const INERTIAL_CAV_THRESHOLD_MI: f64 = 0.3;
    pub const SOUND_SPEED_BLOOD_M_S: f64 = 1540.0;

    #[must_use]
    pub fn from_throat_dp(dp_throat_pa: f64, blood_density: f64) -> Self {
        let v_rms = (2.0 * dp_throat_pa.max(0.0) / blood_density.max(1.0)).sqrt();
        let mi_equiv = v_rms / Self::SOUND_SPEED_BLOOD_M_S;
        Self {
            mi_equiv,
            fda_mi_compliant: mi_equiv < Self::FDA_MI_LIMIT,
            therapeutic_mi_lower: Self::INERTIAL_CAV_THRESHOLD_MI,
            in_therapeutic_window: (Self::INERTIAL_CAV_THRESHOLD_MI..Self::FDA_MI_LIMIT)
                .contains(&mi_equiv),
        }
    }
}

crate::impl_metadata!(FdaCavitationCompliance, "FdaCavitationCompliance");