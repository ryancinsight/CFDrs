//! Physical field identity and display labels.

/// Physical field quantity from a CFD simulation result.
///
/// Each variant identifies a distinct observable mapped onto schematic
/// channels or nodes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnalysisField {
    /// Static pressure `[Pa]`, a node quantity interpolated to edges for display.
    Pressure,
    /// Wall shear stress `[Pa]`, derived from the wall-normal velocity gradient.
    WallShearStress,
    /// Mean cross-sectional velocity `[m/s]`.
    Velocity,
    /// Volumetric flow rate `[m^3/s]`.
    FlowRate,
    /// Apparent dynamic viscosity `[Pa s]`.
    Viscosity,
    /// User-defined field with its display label.
    Custom(String),
}

impl AnalysisField {
    /// Return the axis or legend label including physical units.
    #[must_use]
    pub fn label(&self) -> &str {
        match self {
            Self::Pressure => "Pressure [Pa]",
            Self::WallShearStress => "Wall Shear Stress [Pa]",
            Self::Velocity => "Velocity [m/s]",
            Self::FlowRate => "Flow Rate [m³/s]",
            Self::Viscosity => "Viscosity [Pa·s]",
            Self::Custom(label) => label.as_str(),
        }
    }
}
