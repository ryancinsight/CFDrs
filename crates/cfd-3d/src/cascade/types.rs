use aequitas::systems::si::quantities::{Length, Pressure, Velocity, VolumetricFlowRate};
use serde::{Deserialize, Serialize};

// ── Configuration ─────────────────────────────────────────────────────────────

/// Specification for a single channel in the cascade.
#[derive(Debug, Clone)]
pub struct CascadeChannelSpec {
    /// Identifier matching the `NetworkBlueprint` channel ID.
    pub id: String,
    /// Axial length \[m].
    pub length: Length,
    /// Channel width \[m] (cross-stream, varies for venturi).
    pub width: Length,
    /// Channel height \[m] (constant for rectangular ducts).
    pub height: Length,
    /// Assigned volumetric flow rate [m³/s] from 1D presolver.
    pub flow_rate_m3_s: VolumetricFlowRate,
    /// Whether this channel contains a Venturi throat.
    pub is_venturi_throat: bool,
    /// Throat width \[m] (only when `is_venturi_throat`).
    pub throat_width: Option<Length>,
    /// Local hematocrit [-] from Zweifach-Fung junction routing.
    ///
    /// When set, the solver creates a `CassonBlood::with_hematocrit()` model
    /// for this channel, adjusting yield stress (Chien 1970: τ_y ∝ HCT³) and
    /// infinite-shear viscosity (Quemada 1978) to reflect the local RBC
    /// concentration after upstream cell separation.
    ///
    /// If `None`, the feed fluid is used unchanged.
    pub local_hematocrit: Option<f64>,
}

/// Configuration for the 3D cascade solver.
#[derive(Debug, Clone)]
pub struct CascadeConfig3D {
    /// Reference pressure at outlets \[Pa].
    pub outlet_pressure: Pressure,
    /// Mesh resolution: (axial, cross-width, cross-height).
    pub resolution: (usize, usize, usize),
    /// Maximum Picard iterations for non-Newtonian viscosity coupling.
    pub max_picard_iterations: usize,
    /// Picard convergence tolerance (relative viscosity change).
    pub picard_tolerance: f64,
}

impl Default for CascadeConfig3D {
    fn default() -> Self {
        Self {
            outlet_pressure: Pressure::from_base(0.0),
            resolution: (40, 8, 8),
            max_picard_iterations: 10,
            picard_tolerance: 1e-3,
        }
    }
}

// ── Results ───────────────────────────────────────────────────────────────────

/// Per-channel result from the 3D solve.
#[derive(Debug, Clone)]
pub struct ChannelResult3D {
    /// Channel identifier.
    pub channel_id: String,
    /// Mean wall shear stress \[Pa].
    pub wall_shear_mean_pa: Pressure,
    /// Maximum wall shear stress \[Pa].
    pub wall_shear_max_pa: Pressure,
    /// Pressure drop inlet → outlet \[Pa].
    pub pressure_drop_pa: Pressure,
    /// Maximum velocity magnitude \[m/s].
    pub max_velocity: Velocity,
    /// Whether the solver converged.
    pub converged: bool,
    /// Number of Picard iterations used.
    pub picard_iterations: usize,
    /// Local hematocrit used for this channel's viscosity model.
    /// Equals `CascadeChannelSpec::local_hematocrit` when set, or the feed
    /// fluid's native hematocrit otherwise.
    pub local_hematocrit: f64,
}

/// Aggregate result for the entire cascade.
#[derive(Debug, Clone)]
pub struct CascadeResult3D {
    /// Per-channel results in the order supplied.
    pub channel_results: Vec<ChannelResult3D>,
    /// Sum of per-channel pressure drops \[Pa].
    pub total_pressure_drop_pa: Pressure,
    /// Channel with the highest wall shear stress.
    pub max_shear_channel_id: String,
    /// Highest wall shear stress across all channels \[Pa].
    pub max_shear_pa: Pressure,
}

#[derive(Debug, Serialize, Deserialize)]
struct CascadeChannelSpecRepr {
    id: String,
    length: f64,
    width: f64,
    height: f64,
    flow_rate_m3_s: f64,
    is_venturi_throat: bool,
    throat_width: Option<f64>,
    local_hematocrit: Option<f64>,
}

impl Serialize for CascadeChannelSpec {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        CascadeChannelSpecRepr {
            id: self.id.clone(),
            length: self.length.into_base(),
            width: self.width.into_base(),
            height: self.height.into_base(),
            flow_rate_m3_s: self.flow_rate_m3_s.into_base(),
            is_venturi_throat: self.is_venturi_throat,
            throat_width: self.throat_width.map(Length::into_base),
            local_hematocrit: self.local_hematocrit,
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for CascadeChannelSpec {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let repr = CascadeChannelSpecRepr::deserialize(deserializer)?;
        Ok(Self {
            id: repr.id,
            length: Length::from_base(repr.length),
            width: Length::from_base(repr.width),
            height: Length::from_base(repr.height),
            flow_rate_m3_s: VolumetricFlowRate::from_base(repr.flow_rate_m3_s),
            is_venturi_throat: repr.is_venturi_throat,
            throat_width: repr.throat_width.map(Length::from_base),
            local_hematocrit: repr.local_hematocrit,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct CascadeConfig3DRepr {
    outlet_pressure: f64,
    resolution: (usize, usize, usize),
    max_picard_iterations: usize,
    picard_tolerance: f64,
}

impl Serialize for CascadeConfig3D {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        CascadeConfig3DRepr {
            outlet_pressure: self.outlet_pressure.into_base(),
            resolution: self.resolution,
            max_picard_iterations: self.max_picard_iterations,
            picard_tolerance: self.picard_tolerance,
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for CascadeConfig3D {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let repr = CascadeConfig3DRepr::deserialize(deserializer)?;
        Ok(Self {
            outlet_pressure: Pressure::from_base(repr.outlet_pressure),
            resolution: repr.resolution,
            max_picard_iterations: repr.max_picard_iterations,
            picard_tolerance: repr.picard_tolerance,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct ChannelResult3DRepr {
    channel_id: String,
    wall_shear_mean_pa: f64,
    wall_shear_max_pa: f64,
    pressure_drop_pa: f64,
    max_velocity: f64,
    converged: bool,
    picard_iterations: usize,
    local_hematocrit: f64,
}

impl Serialize for ChannelResult3D {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        ChannelResult3DRepr {
            channel_id: self.channel_id.clone(),
            wall_shear_mean_pa: self.wall_shear_mean_pa.into_base(),
            wall_shear_max_pa: self.wall_shear_max_pa.into_base(),
            pressure_drop_pa: self.pressure_drop_pa.into_base(),
            max_velocity: self.max_velocity.into_base(),
            converged: self.converged,
            picard_iterations: self.picard_iterations,
            local_hematocrit: self.local_hematocrit,
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for ChannelResult3D {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let repr = ChannelResult3DRepr::deserialize(deserializer)?;
        Ok(Self {
            channel_id: repr.channel_id,
            wall_shear_mean_pa: Pressure::from_base(repr.wall_shear_mean_pa),
            wall_shear_max_pa: Pressure::from_base(repr.wall_shear_max_pa),
            pressure_drop_pa: Pressure::from_base(repr.pressure_drop_pa),
            max_velocity: Velocity::from_base(repr.max_velocity),
            converged: repr.converged,
            picard_iterations: repr.picard_iterations,
            local_hematocrit: repr.local_hematocrit,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct CascadeResult3DRepr {
    channel_results: Vec<ChannelResult3D>,
    total_pressure_drop_pa: f64,
    max_shear_channel_id: String,
    max_shear_pa: f64,
}

impl Serialize for CascadeResult3D {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        CascadeResult3DRepr {
            channel_results: self.channel_results.clone(),
            total_pressure_drop_pa: self.total_pressure_drop_pa.into_base(),
            max_shear_channel_id: self.max_shear_channel_id.clone(),
            max_shear_pa: self.max_shear_pa.into_base(),
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for CascadeResult3D {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let repr = CascadeResult3DRepr::deserialize(deserializer)?;
        Ok(Self {
            channel_results: repr.channel_results,
            total_pressure_drop_pa: Pressure::from_base(repr.total_pressure_drop_pa),
            max_shear_channel_id: repr.max_shear_channel_id,
            max_shear_pa: Pressure::from_base(repr.max_shear_pa),
        })
    }
}
