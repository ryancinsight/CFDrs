use super::{EdgeId, NodeId};
use crate::geometry::metadata::{
    ChannelPathMetadata, JunctionGeometryMetadata, NodeLayoutMetadata, VenturiGeometryMetadata,
};
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
            Self::Circular { diameter_m } => std::f64::consts::PI * (diameter_m / 2.0).powi(2),
            Self::Rectangular { width_m, height_m } => width_m * height_m,
        }
    }

    /// Bounding-box dimensions `(width, height)` [m].
    ///
    /// - Circular: `(d, d)`
    /// - Rectangular: `(w, h)`
    #[must_use]
    pub fn dims(&self) -> (f64, f64) {
        match self {
            Self::Circular { diameter_m } => (*diameter_m, *diameter_m),
            Self::Rectangular { width_m, height_m } => (*width_m, *height_m),
        }
    }

    /// Fully-developed Poiseuille wall shear rate [1/s] for a given
    /// mean velocity `u_mean` [m/s].
    ///
    /// ## Theorem — Analytical Wall Shear Rate
    ///
    /// For fully-developed Poiseuille flow in a straight duct:
    ///
    /// - **Rectangular** (Boussinesq 1868): `γ̇_w = 6 u_mean / h`
    ///   where `h` = channel height (smallest dimension governs shear).
    ///
    /// - **Circular** (Hagen-Poiseuille): `γ̇_w = 8 u_mean / d`
    ///   derived from `τ_w = μ · 8V/(πR³·A) = μ · 8u/d`.
    #[must_use]
    pub fn wall_shear_rate(&self, u_mean: f64) -> f64 {
        match self {
            Self::Rectangular { height_m, .. } => 6.0 * u_mean / height_m.max(1e-18),
            Self::Circular { diameter_m } => 8.0 * u_mean / diameter_m.max(1e-18),
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layout: Option<NodeLayoutMetadata>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub junction_geometry: Option<JunctionGeometryMetadata>,
    #[serde(skip)]
    pub metadata: Option<crate::geometry::metadata::MetadataContainer>,
}

impl NodeSpec {
    #[must_use]
    pub fn new(id: impl Into<String>, kind: NodeKind) -> Self {
        Self {
            id: NodeId::new(id),
            kind,
            layout: None,
            junction_geometry: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn with_layout(mut self, layout: NodeLayoutMetadata) -> Self {
        self.layout = Some(layout);
        self
    }

    #[must_use]
    pub fn with_junction_geometry(mut self, geometry: JunctionGeometryMetadata) -> Self {
        self.junction_geometry = Some(geometry);
        self
    }

    #[must_use]
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        meta: T,
    ) -> Self {
        if self.metadata.is_none() {
            self.metadata = Some(crate::geometry::metadata::MetadataContainer::new());
        }
        self.metadata
            .as_mut()
            .expect("metadata container must exist")
            .insert(meta);
        self
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeKind {
    Pipe,
    Valve,
    Pump,
}

/// Channel geometry shape — determines which resistance model the 1D solver
/// dispatches. Defaults to `Straight` for backward compatibility.
///
/// # Theorem — Dean Flow Correction
///
/// In curved channels the secondary (Dean) flow enhances friction:
///
/// `f_curved / f_straight = 1 + 0.033 (log₁₀ De)⁴`  (Mishra & Gupta 1979)
///
/// where `De = Re √(D_h / (2 R_c))` is the Dean number.
///
/// **Proof sketch**: The centripetal acceleration in a curved duct creates a
/// pressure gradient ∂p/∂r ≈ ρu²/R that drives counter-rotating Dean vortices,
/// increasing wall shear stress beyond the Hagen-Poiseuille value.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ChannelShape {
    /// Regular straight duct — Hagen-Poiseuille / Shah-London resistance only.
    Straight,
    /// Serpentine channel with 180° U-turns — triggers Dean flow corrections
    /// and bend minor-loss K-factors in the 1D solver.
    Serpentine {
        /// Number of straight segments between turns.
        segments: usize,
        /// Radius of curvature at U-turn bends \[m\].
        bend_radius_m: f64,
    },
}

impl Default for ChannelShape {
    fn default() -> Self {
        Self::Straight
    }
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
    /// Channel geometry shape — straight vs serpentine w/ Dean corrections.
    #[serde(default)]
    pub channel_shape: ChannelShape,
    pub resistance: f64,
    pub quad_coeff: f64,
    // Component properties
    pub valve_cv: Option<f64>,
    pub pump_max_flow: Option<f64>,
    pub pump_max_pressure: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<ChannelPathMetadata>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub venturi_geometry: Option<VenturiGeometryMetadata>,
    #[serde(skip)]
    pub metadata: Option<crate::geometry::metadata::MetadataContainer>,
}

impl ChannelSpec {
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        meta: T,
    ) -> Self {
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
            channel_shape: ChannelShape::Straight,
            resistance,
            quad_coeff,
            valve_cv: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            path: None,
            venturi_geometry: None,
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
            channel_shape: ChannelShape::Straight,
            resistance,
            quad_coeff,
            valve_cv: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            path: None,
            venturi_geometry: None,
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
            channel_shape: ChannelShape::Straight,
            resistance: 0.0,
            quad_coeff: 0.0,
            valve_cv: Some(cv),
            pump_max_flow: None,
            pump_max_pressure: None,
            path: None,
            venturi_geometry: None,
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
            channel_shape: ChannelShape::Straight,
            resistance: 0.0,
            quad_coeff: 0.0,
            valve_cv: None,
            pump_max_flow: Some(max_flow),
            pump_max_pressure: Some(max_pressure),
            path: None,
            venturi_geometry: None,
            metadata: None,
        }
    }

    #[must_use]
    pub fn with_path(mut self, path: ChannelPathMetadata) -> Self {
        self.path = Some(path);
        self
    }

    #[must_use]
    pub fn with_venturi_geometry(mut self, geometry: VenturiGeometryMetadata) -> Self {
        self.venturi_geometry = Some(geometry);
        self
    }
}
