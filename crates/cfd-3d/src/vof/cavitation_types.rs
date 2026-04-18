//! Cavitation-VOF Configuration and Statistics Types
//!
//! Data types for configuring the cavitation-VOF solver and reporting
//! cavitation statistics.

use crate::vof::config::VofConfig;
use cfd_core::physics::cavitation::{damage::CavitationDamage, models::CavitationModel};

use super::bubble_dynamics::BubbleDynamicsConfig;
use cfd_core::physics::fluid::BloodModel;

/// Cavitation-VOF solver configuration
#[derive(Debug, Clone)]
pub struct CavitationVofConfig {
    /// Base VOF configuration
    pub vof_config: VofConfig,
    /// Cavitation model selection
    pub cavitation_model: CavitationModel<f64>,
    /// Cavitation damage model
    pub damage_model: Option<CavitationDamage<f64>>,
    /// Rayleigh-Plesset bubble dynamics
    pub bubble_dynamics: Option<BubbleDynamicsConfig>,
    /// Passive scalar tracking for cavitation nuclei cascade
    pub nuclei_transport:
        Option<cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig<f64>>,
    /// Cavitation inception threshold
    pub inception_threshold: f64,
    /// Maximum void fraction allowed
    pub max_void_fraction: f64,
    /// Cavitation relaxation time (for numerical stability)
    pub relaxation_time: f64,
    /// Vapor pressure (Pa)
    pub vapor_pressure: f64,
    /// Liquid density (kg/m³)
    pub liquid_density: f64,
    /// Liquid blood rheology model
    pub liquid_blood_model: BloodModel<f64>,
    /// Vapor density (kg/m³)
    pub vapor_density: f64,
    /// Speed of sound in liquid (m/s)
    pub sound_speed: f64,
}

/// Cavitation statistics
#[derive(Debug, Clone)]
pub struct CavitationStatistics {
    /// Fraction of cells experiencing cavitation
    pub cavitation_fraction: f64,
    /// Total void fraction in domain
    pub total_void_fraction: f64,
    /// Maximum void fraction
    pub max_void_fraction: f64,
    /// Maximum accumulated damage
    pub max_damage: f64,
    /// Number of cavitating cells
    pub cavitating_cells: usize,
    /// Total number of cells
    pub total_cells: usize,
}

impl std::fmt::Display for CavitationStatistics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cavitation Statistics:\n\
             Cavitation Fraction: {:.3} ({}/{} cells)\n\
             Total Void Fraction: {:.6}\n\
             Maximum Void Fraction: {:.6}\n\
             Maximum Damage: {:.2e}",
            self.cavitation_fraction,
            self.cavitating_cells,
            self.total_cells,
            self.total_void_fraction,
            self.max_void_fraction,
            self.max_damage
        )
    }
}
