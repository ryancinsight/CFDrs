//! Integration with the scheme library for 2D schematic visualization.
//!
//! This module provides functionality to:
//! - Import microfluidic network designs from 2D schematics
//! - Export 1D simulation results to 2D visualizations
//! - Map between schematic components and simulation elements
//!
//! When the `scheme-integration` feature is disabled, this module provides
//! stub implementations that return appropriate errors.

#[cfg(feature = "scheme-integration")]
use scheme::{
    geometry::{
        types::{ChannelSystem, Channel, Node, Point2D, ChannelType},
        SplitType,
        generator::create_geometry,
    },
    config::{GeometryConfig, ChannelTypeConfig},
    visualizations::schematic::plot_geometry,
    error::SchemeError as SchemeLibError,
};

/// Trait for converting between scheme representations and our network model
pub trait SchemeConversion {
    /// Import a network from a scheme design
    fn from_scheme(scheme_data: &str) -> Result<Self, SchemeError>
    where
        Self: Sized;
    
    /// Export the network to a scheme representation
    fn to_scheme(&self) -> Result<String, SchemeError>;
    
    #[cfg(feature = "scheme-integration")]
    /// Import from a ChannelSystem
    fn from_channel_system(system: &ChannelSystem) -> Result<Self, SchemeError>
    where
        Self: Sized;
    
    #[cfg(feature = "scheme-integration")]
    /// Export to a ChannelSystem
    fn to_channel_system(&self) -> Result<ChannelSystem, SchemeError>;
}

/// Errors that can occur during scheme conversion
#[derive(Debug, thiserror::Error)]
pub enum SchemeError {
    /// Invalid scheme format
    #[error("Invalid scheme format: {0}")]
    InvalidFormat(String),
    
    /// Unsupported component type
    #[error("Unsupported component: {0}")]
    UnsupportedComponent(String),
    
    /// Conversion error
    #[error("Conversion error: {0}")]
    ConversionError(String),
    
    /// Scheme library error
    #[cfg(feature = "scheme-integration")]
    #[error("Scheme library error: {0}")]
    SchemeLibError(#[from] SchemeLibError),
    
    /// JSON serialization error
    #[error("JSON error: {0}")]
    JsonError(#[from] serde_json::Error),

    /// Feature disabled error
    #[error("Feature disabled: {0}")]
    FeatureDisabled(String),
}

/// Component mapping between scheme and simulation
pub struct ComponentMapping {
    /// Scheme component ID
    pub scheme_id: String,
    /// Simulation component ID
    pub sim_id: String,
    /// Component type
    pub component_type: ComponentType,
}

/// Types of components that can be mapped
#[derive(Debug, Clone, PartialEq)]
pub enum ComponentType {
    /// Channel or pipe
    Channel,
    /// Pump
    Pump,
    /// Valve
    Valve,
    /// Junction or mixer (bifurcation/trifurcation)
    Junction(JunctionType),
    /// Sensor
    Sensor,
    /// Inlet/outlet port
    Port,
}

/// Types of junctions based on scheme's split types
#[derive(Debug, Clone, PartialEq)]
pub enum JunctionType {
    /// Two-way split
    Bifurcation,
    /// Three-way split
    Trifurcation,
    /// Custom n-way split
    Custom(usize),
}

#[cfg(feature = "scheme-integration")]
impl From<SplitType> for JunctionType {
    fn from(split: SplitType) -> Self {
        match split {
            SplitType::Bifurcation => JunctionType::Bifurcation,
            SplitType::Trifurcation => JunctionType::Trifurcation,
        }
    }
}

/// Layout information for 2D visualization
pub struct SchematicLayout {
    /// Component positions
    pub positions: Vec<(String, f64, f64)>,
    /// Connection paths
    pub connections: Vec<ConnectionPath>,
    /// Layout bounds
    pub bounds: (f64, f64, f64, f64), // min_x, min_y, max_x, max_y
}

/// Path information for a connection
pub struct ConnectionPath {
    /// Source component ID
    pub source: String,
    /// Target component ID
    pub target: String,
    /// Path points
    pub points: Vec<(f64, f64)>,
    /// Channel type (straight, serpentine, etc.)
    pub channel_type: ChannelPathType,
}

/// Channel path types based on scheme's `ChannelType`
#[derive(Debug, Clone)]
pub enum ChannelPathType {
    /// Straight line channel
    Straight,
    /// Serpentine (S-shaped) channel
    Serpentine,
    /// Arc channel
    Arc,
    /// Tapered channel
    Frustum,
}

#[cfg(feature = "scheme-integration")]
impl From<&ChannelType> for ChannelPathType {
    fn from(channel_type: &ChannelType) -> Self {
        match channel_type {
            ChannelType::Straight => ChannelPathType::Straight,
            ChannelType::SmoothStraight { .. } => ChannelPathType::Straight,
            ChannelType::Serpentine { .. } => ChannelPathType::Serpentine,
            ChannelType::Arc { .. } => ChannelPathType::Arc,
            ChannelType::Frustum { .. } => ChannelPathType::Frustum,
        }
    }
}

/// Helper functions for scheme integration
#[cfg(feature = "scheme-integration")]
pub mod helpers {
    use super::*;
    
    /// Create a simple bifurcation network schematic
    pub fn create_bifurcation_schematic(
        width: f64,
        height: f64,
        levels: usize,
    ) -> Result<ChannelSystem, SchemeError> {
        let config = GeometryConfig::default();
        let splits = vec![SplitType::Bifurcation; levels];
        
        Ok(create_geometry(
            (width, height),
            &splits,
            &config,
            &ChannelTypeConfig::AllStraight,
        ))
    }
    
    /// Export a channel system to PNG
    pub fn export_to_png(
        system: &ChannelSystem,
        filename: &str,
    ) -> Result<(), SchemeError> {
        plot_geometry(system, filename)?;
        Ok(())
    }
    
    /// Convert scheme nodes to network nodes
    pub fn convert_nodes(scheme_nodes: &[Node]) -> Vec<(usize, Point2D)> {
        scheme_nodes
            .iter()
            .map(|node| (node.id.clone(), node.point))
            .collect()
    }
    
    /// Extract channel paths
    pub fn extract_channel_paths(channels: &[Channel]) -> Vec<ConnectionPath> {
        channels
            .iter()
            .map(|channel| {
                let path_points = match &channel.channel_type {
                    ChannelType::Straight => {
                        vec![channel.start, channel.end]
                    }
                    ChannelType::SmoothStraight { path } |
                    ChannelType::Serpentine { path } |
                    ChannelType::Arc { path } => path,
                    ChannelType::Frustum { path, .. } => path,
                };
                
                ConnectionPath {
                    source: format!("node_{}", channel.start_node),
                    target: format!("node_{}", channel.end_node),
                    points: path_points,
                    channel_type: (&channel.channel_type).into(),
                }
            })
            .collect()
    }
}

/// Fallback implementations when scheme integration is disabled
#[cfg(not(feature = "scheme-integration"))]
pub mod helpers {
    use super::SchemeError;

    /// Create a simple bifurcation network schematic (fallback)
    pub fn create_bifurcation_schematic(
        _width: f64,
        _height: f64,
        _levels: usize,
    ) -> Result<(), SchemeError> {
        Err(SchemeError::FeatureDisabled(
            "scheme-integration feature is not enabled".to_string()
        ))
    }

    /// Export a channel system to PNG (fallback)
    pub fn export_to_png(
        _system: &(),
        _filename: &str,
    ) -> Result<(), SchemeError> {
        Err(SchemeError::FeatureDisabled(
            "scheme-integration feature is not enabled".to_string()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_component_type() {
        let comp_type = ComponentType::Channel;
        assert_eq!(comp_type, ComponentType::Channel);
    }
    
    #[test]
    fn test_junction_type() {
        let junction = JunctionType::Bifurcation;
        assert_eq!(junction, JunctionType::Bifurcation);
    }
    
    #[cfg(feature = "scheme-integration")]
    #[test]
    fn test_split_type_conversion() {
        assert_eq!(
            JunctionType::from(SplitType::Bifurcation),
            JunctionType::Bifurcation
        );
        assert_eq!(
            JunctionType::from(SplitType::Trifurcation),
            JunctionType::Trifurcation
        );
    }
}