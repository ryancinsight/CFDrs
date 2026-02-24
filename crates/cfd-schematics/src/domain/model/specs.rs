use super::{EdgeId, NodeId};
use serde::{Deserialize, Serialize};

/// Cross-section geometry specification for a channel
///
/// Carries the physical geometry needed by the 1D solver to compute
/// hydraulic resistance. This lives in `cfd-schematics` so that examples
/// can fully specify a network via `ChannelSpec` without importing
/// `cfd_1d::channel` types.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum CrossSectionSpec {
    /// Circular cross-section (e.g. tubing, blood vessels)
    Circular {
        /// Inner diameter [m]
        diameter_m: f64,
    },
    /// Rectangular cross-section (e.g. PDMS microfluidic channels)
    Rectangular {
        /// Width [m]
        width_m: f64,
        /// Height [m]
        height_m: f64,
    },
}

impl CrossSectionSpec {
    /// Hydraulic diameter [m]: `D_h = 4A / P`
    ///
    /// - Circular: `D_h = d`
    /// - Rectangular: `D_h = 2wh / (w + h)`
    #[must_use]
    pub fn hydraulic_diameter(&self) -> f64 {
        match self {
            Self::Circular { diameter_m } => *diameter_m,
            Self::Rectangular { width_m, height_m } => {
                2.0 * width_m * height_m / (width_m + height_m)
            }
        }
    }

    /// Cross-sectional area [m²]
    #[must_use]
    pub fn area(&self) -> f64 {
        match self {
            Self::Circular { diameter_m } => {
                std::f64::consts::PI * (diameter_m / 2.0).powi(2)
            }
            Self::Rectangular { width_m, height_m } => width_m * height_m,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeKind {
    Inlet,
    Outlet,
    Reservoir,
    Junction,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeSpec {
    pub id: NodeId,
    pub kind: NodeKind,
}

impl NodeSpec {
    #[must_use]
    pub fn new(id: impl Into<String>, kind: NodeKind) -> Self {
        Self {
            id: NodeId::new(id),
            kind,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeKind {
    Pipe,
    Valve,
    Pump,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelSpec {
    pub id: EdgeId,
    pub kind: EdgeKind,
    pub from: NodeId,
    pub to: NodeId,
    pub length_m: f64,
    /// Physical cross-section geometry — used by `cfd-1d` to compute hydraulic resistance.
    /// Replaces the old `diameter_m` field with a typed, extensible spec.
    pub cross_section: CrossSectionSpec,
    pub resistance: f64,
    pub quad_coeff: f64,
    // Component properties
    pub valve_cv: Option<f64>,
    pub pump_max_flow: Option<f64>,
    pub pump_max_pressure: Option<f64>,
    #[serde(skip)]
    pub metadata: Option<crate::geometry::metadata::MetadataContainer>,
}

impl ChannelSpec {
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(mut self, meta: T) -> Self {
        if self.metadata.is_none() {
            self.metadata = Some(crate::geometry::metadata::MetadataContainer::new());
        }
        self.metadata.as_mut().unwrap().insert(meta);
        self
    }

    /// Create a circular-cross-section pipe channel spec.
    #[must_use]
    pub fn new_pipe(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        length_m: f64,
        diameter_m: f64,
        resistance: f64,
        quad_coeff: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Pipe,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m,
            cross_section: CrossSectionSpec::Circular { diameter_m },
            resistance,
            quad_coeff,
            valve_cv: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            metadata: None,
        }
    }

    /// Create a rectangular-cross-section pipe channel spec.
    #[must_use]
    pub fn new_pipe_rect(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        length_m: f64,
        width_m: f64,
        height_m: f64,
        resistance: f64,
        quad_coeff: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Pipe,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m,
            cross_section: CrossSectionSpec::Rectangular { width_m, height_m },
            resistance,
            quad_coeff,
            valve_cv: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn new_valve(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        cv: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Valve,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m: 0.0,
            cross_section: CrossSectionSpec::Circular { diameter_m: 0.0 },
            resistance: 0.0,
            quad_coeff: 0.0,
            valve_cv: Some(cv),
            pump_max_flow: None,
            pump_max_pressure: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn new_pump(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        max_flow: f64,
        max_pressure: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Pump,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m: 0.0,
            cross_section: CrossSectionSpec::Circular { diameter_m: 0.0 },
            resistance: 0.0,
            quad_coeff: 0.0,
            valve_cv: None,
            pump_max_flow: Some(max_flow),
            pump_max_pressure: Some(max_pressure),
            metadata: None,
        }
    }
}
