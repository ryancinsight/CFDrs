use std::collections::{HashMap, HashSet, VecDeque};

use cfd_core::error::{Error, Result};
use cfd_schematics::domain::model::{NetworkBlueprint, NodeKind};
use cfd_schematics::geometry::metadata::{ChannelVenturiSpec, GeometryAuthoringProvenance};

const ENDPOINT_TOLERANCE_MM: f64 = 1.0e-3;
const LENGTH_REL_TOL: f64 = 5.0e-2;
const LENGTH_ABS_TOL_M: f64 = 1.0e-6;

/// Validate that a blueprint is admissible for the 1D network solve and that
/// its explicit figure geometry is consistent with the solver graph contract.
///
/// # Errors
/// Returns an [`Error::InvalidConfiguration`] describing the first violated
/// structural or geometry-to-graph invariant.
pub fn validate_blueprint_for_1d_solve(blueprint: &NetworkBlueprint) -> Result<()> {
    cfd_schematics::domain::rules::BlueprintValidator::validate(blueprint)?;

    if blueprint.inlet_count() == 0 || blueprint.outlet_count() == 0 {
        return Err(Error::InvalidConfiguration(
            "Blueprint requires at least one inlet and one outlet".to_string(),
        ));
    }

    let node_points: HashMap<&str, (f64, f64)> = blueprint
        .nodes
        .iter()
        .map(|node| (node.id.as_str(), node.point))
        .collect();

    let mut adjacency: HashMap<&str, Vec<&str>> = HashMap::with_capacity(blueprint.nodes.len());
    for channel in &blueprint.channels {
        adjacency
            .entry(channel.from.as_str())
            .or_default()
            .push(channel.to.as_str());

        let from = node_points
            .get(channel.from.as_str())
            .copied()
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "Channel '{}' references unknown source node '{}'",
                    channel.id.as_str(),
                    channel.from.as_str()
                ))
            })?;
        let to = node_points
            .get(channel.to.as_str())
            .copied()
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "Channel '{}' references unknown target node '{}'",
                    channel.id.as_str(),
                    channel.to.as_str()
                ))
            })?;

        if channel.path.len() < 2 {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' is missing an explicit routed path",
                channel.id.as_str()
            )));
        }

        let start = channel.path[0];
        let end = *channel.path.last().expect("path len >= 2");
        if point_distance_mm(start, from) > ENDPOINT_TOLERANCE_MM {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' path start does not match node '{}'",
                channel.id.as_str(),
                channel.from.as_str()
            )));
        }
        if point_distance_mm(end, to) > ENDPOINT_TOLERANCE_MM {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' path end does not match node '{}'",
                channel.id.as_str(),
                channel.to.as_str()
            )));
        }

        let polyline_length_m = polyline_length_mm(&channel.path) * 1.0e-3;
        let length_tolerance = LENGTH_ABS_TOL_M.max(channel.length_m.abs() * LENGTH_REL_TOL);

        // Skip tight length checks for symbolic schematics (e.g. from selective_wrapper)
        let is_symbolic = blueprint
            .metadata
            .as_ref()
            .and_then(|meta| meta.get::<GeometryAuthoringProvenance>())
            .is_some_and(|prov| prov.source.contains("selective_wrapper"));

        if !is_symbolic && (polyline_length_m - channel.length_m).abs() > length_tolerance {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' path length {:.6e} m disagrees with declared length {:.6e} m",
                channel.id.as_str(),
                polyline_length_m,
                channel.length_m
            )));
        }

        let serial_venturi_throats = channel
            .metadata
            .as_ref()
            .and_then(|meta| meta.get::<ChannelVenturiSpec>())
            .map_or(0, |spec| spec.n_throats);
        if serial_venturi_throats > 0 && channel.venturi_geometry.is_none() {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' declares venturi throats but has no explicit venturi geometry",
                channel.id.as_str()
            )));
        }
    }

    validate_connected_inlet_to_outlet(blueprint, &adjacency)?;
    validate_duplicate_channel_geometry(blueprint)?;

    Ok(())
}

fn validate_connected_inlet_to_outlet(
    blueprint: &NetworkBlueprint,
    adjacency: &HashMap<&str, Vec<&str>>,
) -> Result<()> {
    let inlet_ids: Vec<&str> = blueprint
        .nodes
        .iter()
        .filter(|node| node.kind == NodeKind::Inlet)
        .map(|node| node.id.as_str())
        .collect();
    let outlet_ids: HashSet<&str> = blueprint
        .nodes
        .iter()
        .filter(|node| node.kind == NodeKind::Outlet)
        .map(|node| node.id.as_str())
        .collect();

    let mut queue: VecDeque<&str> = inlet_ids.iter().copied().collect();
    let mut visited = HashSet::new();
    while let Some(node_id) = queue.pop_front() {
        if !visited.insert(node_id) {
            continue;
        }
        if outlet_ids.contains(node_id) {
            return Ok(());
        }
        if let Some(neighbors) = adjacency.get(node_id) {
            for neighbor in neighbors {
                queue.push_back(neighbor);
            }
        }
    }

    Err(Error::InvalidConfiguration(
        "Blueprint has no directed inlet-to-outlet path".to_string(),
    ))
}

fn validate_duplicate_channel_geometry(blueprint: &NetworkBlueprint) -> Result<()> {
    let mut seen = HashSet::with_capacity(blueprint.channels.len());
    for channel in &blueprint.channels {
        let signature = format!(
            "{}|{}|{}",
            channel.from.as_str(),
            channel.to.as_str(),
            rounded_path_signature(&channel.path)
        );
        if !seen.insert(signature) {
            return Err(Error::InvalidConfiguration(format!(
                "Blueprint contains duplicate overlapping channel geometry for '{}'",
                channel.id.as_str()
            )));
        }
    }
    Ok(())
}

fn point_distance_mm(left: (f64, f64), right: (f64, f64)) -> f64 {
    (left.0 - right.0).hypot(left.1 - right.1)
}

fn polyline_length_mm(path: &[(f64, f64)]) -> f64 {
    path.windows(2)
        .map(|segment| (segment[1].0 - segment[0].0).hypot(segment[1].1 - segment[0].1))
        .sum()
}

fn rounded_path_signature(path: &[(f64, f64)]) -> String {
    path.iter()
        .map(|(x, y)| format!("{x:.6}:{y:.6}"))
        .collect::<Vec<_>>()
        .join(";")
}
