//! Hydraulic resistance models for 1D CFD networks.
//!
//! This module provides comprehensive resistance modeling for various
//! microfluidic components and flow conditions, including analytical
//! solutions and empirical correlations.

pub mod calculator;
pub mod factory;
pub mod geometry;
pub mod models;

// Re-export main types
pub use calculator::ResistanceCalculator;
pub use factory::ResistanceModelFactory;
pub use geometry::ChannelGeometry;
pub use models::CombinationMethod;
pub use models::{
    BendType, DarcyWeisbachModel, EntranceEffectsModel, ExpansionType, FlowConditions,
    HagenPoiseuilleModel, RectangularChannelModel, ResistanceModel, SerpentineAnalysis,
    SerpentineCrossSection, SerpentineModel, VenturiAnalysis, VenturiGeometry, VenturiModel,
};
