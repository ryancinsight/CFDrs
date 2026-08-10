use crate::domain::model::{NetworkBlueprint, NodeKind};
use cfd_core::error::{Error, Result};
use std::collections::HashSet;

/// Stateless validator enforcing structural integrity of a network blueprint.
pub struct BlueprintValidator;

impl BlueprintValidator {
    /// Validate blueprint structure: non-empty node/channel sets, at least one
    /// inlet and outlet, unique identifiers, and resolvable, physically
    /// positive channel geometry.
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidConfiguration` on the first violated rule.
    pub fn validate(blueprint: &NetworkBlueprint) -> Result<()> {
        if blueprint.nodes.is_empty() {
            return Err(Error::InvalidConfiguration(
                "Blueprint contains no nodes".to_string(),
            ));
        }

        if blueprint.channels.is_empty() {
            return Err(Error::InvalidConfiguration(
                "Blueprint contains no channels".to_string(),
            ));
        }

        let inlet_count = blueprint
            .nodes
            .iter()
            .filter(|node| node.kind == NodeKind::Inlet)
            .count();
        let outlet_count = blueprint
            .nodes
            .iter()
            .filter(|node| node.kind == NodeKind::Outlet)
            .count();

        if inlet_count == 0 {
            return Err(Error::InvalidConfiguration(
                "Blueprint requires at least one inlet".to_string(),
            ));
        }

        if outlet_count == 0 {
            return Err(Error::InvalidConfiguration(
                "Blueprint requires at least one outlet".to_string(),
            ));
        }

        let mut node_ids = HashSet::new();
        for node in &blueprint.nodes {
            if !node_ids.insert(node.id.as_str().to_string()) {
                return Err(Error::InvalidConfiguration(format!(
                    "Duplicate node id: {}",
                    node.id.as_str()
                )));
            }
        }

        let mut edge_ids = HashSet::new();
        for channel in &blueprint.channels {
            if !edge_ids.insert(channel.id.as_str().to_string()) {
                return Err(Error::InvalidConfiguration(format!(
                    "Duplicate channel id: {}",
                    channel.id.as_str()
                )));
            }

            if !node_ids.contains(channel.from.as_str()) {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' references unknown source node '{}'",
                    channel.id.as_str(),
                    channel.from.as_str()
                )));
            }

            if !node_ids.contains(channel.to.as_str()) {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' references unknown target node '{}'",
                    channel.id.as_str(),
                    channel.to.as_str()
                )));
            }

            let length_m = channel.length_m.into_base();
            if length_m <= 0.0 {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' has non-positive length: {}",
                    channel.id.as_str(),
                    length_m
                )));
            }

            let hydraulic_diameter_m = channel.cross_section.hydraulic_diameter().into_base();
            if hydraulic_diameter_m <= 0.0 {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' has non-positive hydraulic diameter: {}",
                    channel.id.as_str(),
                    hydraulic_diameter_m
                )));
            }

            let resistance = channel.resistance.into_base();
            if resistance < 0.0 {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' has negative resistance: {}",
                    channel.id.as_str(),
                    resistance
                )));
            }

            let quad_coeff = channel.quad_coeff.into_base();
            if quad_coeff < 0.0 {
                return Err(Error::InvalidConfiguration(format!(
                    "Channel '{}' has negative quadratic loss coefficient: {}",
                    channel.id.as_str(),
                    quad_coeff
                )));
            }
        }

        Ok(())
    }
}
