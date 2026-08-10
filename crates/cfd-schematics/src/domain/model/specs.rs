use super::{EdgeId, NodeId};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{
    JunctionGeometryMetadata, NodeLayoutMetadata, VenturiGeometryMetadata,
};
use aequitas::systems::si::quantities::{
    Area, HydraulicResistance, Length, Pressure, QuadraticHydraulicResistance, VolumetricFlowRate,
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
        /// Inner diameter as an Eunomia-backed SI length.
        diameter_m: Length<f64>,
    },
    /// Rectangular cross-section (e.g. PDMS microfluidic channels)
    Rectangular {
        /// Width as an Eunomia-backed SI length.
        width_m: Length<f64>,
        /// Height as an Eunomia-backed SI length.
        height_m: Length<f64>,
    },
}

impl CrossSectionSpec {
    /// Hydraulic diameter \[m]: `D_h = 4A / P`
    ///
    /// - Circular: `D_h = d`
    /// - Rectangular: `D_h = 2wh / (w + h)`
    #[must_use]
    pub fn hydraulic_diameter(&self) -> Length<f64> {
        match self {
            Self::Circular { diameter_m } => *diameter_m,
            Self::Rectangular { width_m, height_m } => {
                let width_m = width_m.into_base();
                let height_m = height_m.into_base();
                Length::from_base(2.0 * width_m * height_m / (width_m + height_m))
            }
        }
    }

    /// Cross-sectional area \[m²]
    #[must_use]
    pub fn area(&self) -> Area<f64> {
        match self {
            Self::Circular { diameter_m } => {
                Area::from_base(std::f64::consts::PI * (diameter_m.into_base() / 2.0).powi(2))
            }
            Self::Rectangular { width_m, height_m } => {
                Area::from_base(width_m.into_base() * height_m.into_base())
            }
        }
    }

    /// Bounding-box dimensions `(width, height)` \[m].
    ///
    /// - Circular: `(d, d)`
    /// - Rectangular: `(w, h)`
    #[must_use]
    pub fn dims(&self) -> (Length<f64>, Length<f64>) {
        match self {
            Self::Circular { diameter_m } => (*diameter_m, *diameter_m),
            Self::Rectangular { width_m, height_m } => (*width_m, *height_m),
        }
    }

    /// Fully-developed Poiseuille wall shear rate [1/s] for a given
    /// mean velocity `u_mean` \[m/s].
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
            Self::Rectangular { height_m, .. } => 6.0 * u_mean / height_m.into_base().max(1e-18),
            Self::Circular { diameter_m } => 8.0 * u_mean / diameter_m.into_base().max(1e-18),
        }
    }
}

/// Role of a node in the network topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeKind {
    /// Fluid entry point into the network.
    Inlet,
    /// Fluid exit point from the network.
    Outlet,
    /// Pressure reservoir boundary.
    Reservoir,
    /// Internal branch point connecting multiple edges.
    Junction,
}

/// Node placement and role specification for a network blueprint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeSpec {
    /// Stable identifier of the node.
    pub id: NodeId,
    /// Role of the node in the topology.
    pub kind: NodeKind,
    /// Placement `(x, y)` in schematic coordinates.
    pub point: (f64, f64),
    /// Optional layout metadata for the node.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layout: Option<NodeLayoutMetadata>,
    /// Optional junction geometry metadata for the node.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub junction_geometry: Option<JunctionGeometryMetadata>,
    /// Extension metadata attached to the node.
    #[serde(skip)]
    pub metadata: Option<crate::geometry::metadata::MetadataContainer>,
}

impl NodeSpec {
    /// Construct a node at the origin with the given id and kind.
    #[must_use]
    pub fn new(id: impl Into<String>, kind: NodeKind) -> Self {
        Self::new_at(id, kind, (0.0, 0.0))
    }

    /// Construct a node at the given point with the given id and kind.
    #[must_use]
    pub fn new_at(id: impl Into<String>, kind: NodeKind, point: (f64, f64)) -> Self {
        Self {
            id: NodeId::new(id),
            kind,
            point,
            layout: None,
            junction_geometry: None,
            metadata: None,
        }
    }

    /// Set the placement point.
    #[must_use]
    pub fn with_point(mut self, point: (f64, f64)) -> Self {
        self.point = point;
        self
    }

    /// Set the layout metadata and align the placement point to it.
    #[must_use]
    pub fn with_layout(mut self, layout: NodeLayoutMetadata) -> Self {
        self.point = (layout.x_mm, layout.y_mm);
        self.layout = Some(layout);
        self
    }

    /// Attach junction geometry metadata.
    #[must_use]
    pub fn with_junction_geometry(mut self, geometry: JunctionGeometryMetadata) -> Self {
        self.junction_geometry = Some(geometry);
        self
    }

    /// Attach a typed metadata entry to the node.
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

/// Role of an edge in the network topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeKind {
    /// Passive channel with hydraulic resistance.
    Pipe,
    /// Flow-control valve.
    Valve,
    /// Active pump with flow and pressure limits.
    Pump,
}

/// Channel geometry shape — determines which resistance model the 1D solver
/// dispatches. Defaults to `Straight`.
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
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum ChannelShape {
    /// Regular straight duct — Hagen-Poiseuille / Shah-London resistance only.
    #[default]
    Straight,
    /// Serpentine channel with 180° U-turns - triggers Dean flow corrections
    /// and bend minor-loss K-factors in the 1D solver.
    Serpentine {
        /// Number of straight segments between turns.
        segments: usize,
        /// Radius of curvature at U-turn bends.
        bend_radius_m: Length<f64>,
        /// Waveform type (sine, square, triangular).  Controls the bend
        /// K-factor model: square waves use `BendType::Sharp`, sine and
        /// triangular use `BendType::Smooth` with the specified R/D_h.
        #[serde(default)]
        wave_type: crate::topology::SerpentineWaveType,
    },
}

/// Edge specification connecting two nodes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelSpec {
    /// Stable identifier of the edge.
    pub id: EdgeId,
    /// Role of the edge.
    pub kind: EdgeKind,
    /// Origin node identifier.
    pub from: NodeId,
    /// Destination node identifier.
    pub to: NodeId,
    /// Authored centerline length as an Eunomia-backed SI length.
    pub length_m: Length<f64>,
    /// Physical cross-section geometry — used by `cfd-1d` to compute hydraulic resistance.
    /// Replaces the old `diameter_m` field with a typed, extensible spec.
    pub cross_section: CrossSectionSpec,
    /// Channel geometry shape — straight vs serpentine w/ Dean corrections.
    #[serde(default)]
    pub channel_shape: ChannelShape,
    /// Linear pressure-loss coefficient as an Eunomia-backed SI quantity.
    pub resistance: HydraulicResistance<f64>,
    /// Quadratic pressure-loss coefficient as an Eunomia-backed SI quantity.
    pub quad_coeff: QuadraticHydraulicResistance<f64>,
    // Component properties
    /// Valve quadratic loss coefficient derived from the conventional `Cv` input.
    ///
    /// Eunomia's SI dimensions use integer exponents, so the equivalent
    /// `1 / Cv²` coefficient is the dimensionally exact provider quantity at
    /// this boundary; the square-root pressure form of `Cv` is not representable
    /// as a standalone SI dimension.
    pub valve_loss_coefficient: Option<QuadraticHydraulicResistance<f64>>,
    /// Free-delivery pump flow limit.
    pub pump_max_flow: Option<VolumetricFlowRate<f64>>,
    /// Stall-point pump pressure limit.
    pub pump_max_pressure: Option<Pressure<f64>>,
    /// Visual role used when rendering the edge.
    pub visual_role: Option<crate::geometry::metadata::ChannelVisualRole>,
    /// Authored centerline waypoints in schematic coordinates.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub path: Vec<(f64, f64)>,
    /// Therapy zone the edge belongs to, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub therapy_zone: Option<TherapyZone>,
    /// Venturi geometry metadata for therapy channels.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub venturi_geometry: Option<VenturiGeometryMetadata>,
    /// Extension metadata attached to the edge.
    #[serde(skip)]
    pub metadata: Option<crate::geometry::metadata::MetadataContainer>,
}

impl ChannelSpec {
    /// Attach a typed metadata entry to the channel.
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        meta: T,
    ) -> Self {
        if self.metadata.is_none() {
            self.metadata = Some(crate::geometry::metadata::MetadataContainer::new());
        }
        self.metadata
            .as_mut()
            .expect("structural invariant")
            .insert(meta);
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
            length_m: Length::from_base(length_m),
            cross_section: CrossSectionSpec::Circular {
                diameter_m: Length::from_base(diameter_m),
            },
            channel_shape: ChannelShape::Straight,
            resistance: HydraulicResistance::from_base(resistance),
            quad_coeff: QuadraticHydraulicResistance::from_base(quad_coeff),
            valve_loss_coefficient: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            visual_role: None,
            path: Vec::new(),
            therapy_zone: None,
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
            length_m: Length::from_base(length_m),
            cross_section: CrossSectionSpec::Rectangular {
                width_m: Length::from_base(width_m),
                height_m: Length::from_base(height_m),
            },
            channel_shape: ChannelShape::Straight,
            resistance: HydraulicResistance::from_base(resistance),
            quad_coeff: QuadraticHydraulicResistance::from_base(quad_coeff),
            valve_loss_coefficient: None,
            pump_max_flow: None,
            pump_max_pressure: None,
            visual_role: None,
            path: Vec::new(),
            therapy_zone: None,
            venturi_geometry: None,
            metadata: None,
        }
    }

    /// Create a valve channel spec from a conventional `Cv` loss coefficient.
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
            length_m: Length::from_base(0.0),
            cross_section: CrossSectionSpec::Circular {
                diameter_m: Length::from_base(0.0),
            },
            channel_shape: ChannelShape::Straight,
            resistance: HydraulicResistance::from_base(0.0),
            quad_coeff: QuadraticHydraulicResistance::from_base(0.0),
            valve_loss_coefficient: (cv > 0.0)
                .then(|| QuadraticHydraulicResistance::from_base(1.0 / (cv * cv))),
            pump_max_flow: None,
            pump_max_pressure: None,
            visual_role: None,
            path: Vec::new(),
            therapy_zone: None,
            venturi_geometry: None,
            metadata: None,
        }
    }

    /// Create a pump channel spec with free-delivery flow and stall pressure
    /// limits.
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
            length_m: Length::from_base(0.0),
            cross_section: CrossSectionSpec::Circular {
                diameter_m: Length::from_base(0.0),
            },
            channel_shape: ChannelShape::Straight,
            resistance: HydraulicResistance::from_base(0.0),
            quad_coeff: QuadraticHydraulicResistance::from_base(0.0),
            valve_loss_coefficient: None,
            pump_max_flow: Some(VolumetricFlowRate::from_base(max_flow)),
            pump_max_pressure: Some(Pressure::from_base(max_pressure)),
            visual_role: None,
            path: Vec::new(),
            therapy_zone: None,
            venturi_geometry: None,
            metadata: None,
        }
    }

    /// Set the authored centerline waypoints.
    #[must_use]
    pub fn with_path(mut self, path: Vec<(f64, f64)>) -> Self {
        self.path = path;
        self
    }

    /// Set the visual role used when rendering the edge.
    #[must_use]
    pub fn with_visual_role(mut self, role: crate::geometry::metadata::ChannelVisualRole) -> Self {
        self.visual_role = Some(role);
        self
    }

    /// Attach venturi geometry metadata for therapy channels.
    #[must_use]
    pub fn with_venturi_geometry(mut self, geometry: VenturiGeometryMetadata) -> Self {
        self.venturi_geometry = Some(geometry);
        self
    }

    /// Extract the effective hydraulic width of the channel in meters.
    #[must_use]
    pub fn effective_width_m(&self) -> Length<f64> {
        match self.cross_section {
            CrossSectionSpec::Rectangular { width_m, .. } => width_m,
            CrossSectionSpec::Circular { diameter_m } => diameter_m,
        }
    }

    /// Extract the effective hydraulic height of the channel in meters.
    #[must_use]
    pub fn effective_height_m(&self) -> Length<f64> {
        match self.cross_section {
            CrossSectionSpec::Rectangular { height_m, .. } => height_m,
            CrossSectionSpec::Circular { diameter_m } => diameter_m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::ChannelSpec;

    #[test]
    fn hydraulic_seed_metrics_are_provider_typed() {
        let channel = ChannelSpec::new_pipe("pipe", "in", "out", 0.25, 1.0e-3, 12.0, 3.5);

        assert_eq!(channel.resistance.into_base().to_bits(), 12.0_f64.to_bits());
        assert_eq!(channel.quad_coeff.into_base().to_bits(), 3.5_f64.to_bits());
    }

    #[test]
    fn valve_cv_is_materialized_as_quadratic_loss() {
        let valve = ChannelSpec::new_valve("valve", "in", "out", 0.5);

        assert_eq!(
            valve
                .valve_loss_coefficient
                .expect("positive Cv produces a loss coefficient")
                .into_base()
                .to_bits(),
            4.0_f64.to_bits()
        );
    }

    #[test]
    fn pump_limits_are_provider_typed() {
        let pump = ChannelSpec::new_pump("pump", "in", "out", 2.5e-6, 12_000.0);

        assert_eq!(
            pump.pump_max_flow
                .expect("pump flow limit is present")
                .into_base()
                .to_bits(),
            2.5e-6_f64.to_bits()
        );
        assert_eq!(
            pump.pump_max_pressure
                .expect("pump pressure limit is present")
                .into_base()
                .to_bits(),
            12_000.0_f64.to_bits()
        );
    }
}
