//! Core resistance models for 1D flow calculations.

mod conductance;
mod darcy_weisbach;
mod entrance;
mod hagen_poiseuille;
mod junction_loss;
mod membrane;
mod rectangular;
mod serpentine;
pub(crate) mod traits;
mod venturi;

pub use conductance::{cascade_treatment_flow_fractions, parallel_channel_flow_fractions};
pub use darcy_weisbach::DarcyWeisbachModel;
pub use entrance::{CombinationMethod, EntranceEffectsModel};
pub use hagen_poiseuille::HagenPoiseuilleModel;
pub use junction_loss::{JunctionFlowDirection, JunctionLossModel, JunctionType};
pub use membrane::MembranePoreModel;
pub use rectangular::RectangularChannelModel;
pub use serpentine::{BendType, SerpentineAnalysis, SerpentineCrossSection, SerpentineModel};
pub use traits::{FlowConditions, ResistanceModel};
pub use venturi::{ExpansionType, VenturiAnalysis, VenturiGeometry, VenturiModel};
