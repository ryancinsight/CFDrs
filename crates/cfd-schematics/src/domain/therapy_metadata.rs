use crate::geometry::metadata::Metadata;
use serde::{Deserialize, Serialize};
use std::any::Any;

/// Identifies the routing purpose of a specific channel in sonodynamic therapy.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TherapyZone {
    /// Target zone intended for diseased cells (high cavitation energy).
    CancerTarget,
    /// Bypass zone intended for separating healthy cells (low shear).
    HealthyBypass,
    /// Undifferentiated or pre-separation region.
    MixedFlow,
}

/// Metadata struct for tagging nodes/channels with their therapy routing intent.
/// Implements `cfd_schematics::geometry::metadata::Metadata`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TherapyZoneMetadata {
    pub zone: TherapyZone,
}

impl TherapyZoneMetadata {
    pub fn new(zone: TherapyZone) -> Self {
        Self { zone }
    }
}

impl Metadata for TherapyZoneMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "TherapyZoneMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}
