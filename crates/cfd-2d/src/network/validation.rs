use cfd_core::error::{Error, Result as CfdResult};
use cfd_schematics::domain::model::NetworkBlueprint;

/// Validate that a blueprint is admissible for 2D channel-by-channel
/// projection and meshing.
///
/// # Errors
/// Returns an error when the blueprint fails 1D graph validation or when the
/// routed geometry is unsuitable for 2D projection.
pub fn validate_blueprint_for_2d_projection(blueprint: &NetworkBlueprint) -> CfdResult<()> {
    cfd_1d::validate_blueprint_for_1d_solve(blueprint)?;

    if blueprint.has_unresolved_channel_overlaps() {
        return Err(Error::InvalidConfiguration(format!(
            "Blueprint '{}' contains unresolved routed channel crossings",
            blueprint.name
        )));
    }

    for channel in &blueprint.channels {
        if channel.path.len() < 2 {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' has no explicit routed path for 2D projection",
                channel.id.as_str()
            )));
        }

        let distinct_points = channel
            .path
            .windows(2)
            .filter(|segment| {
                (segment[1].0 - segment[0].0).abs() > f64::EPSILON
                    || (segment[1].1 - segment[0].1).abs() > f64::EPSILON
            })
            .count();
        if distinct_points == 0 {
            return Err(Error::InvalidConfiguration(format!(
                "Channel '{}' has degenerate routed geometry for 2D projection",
                channel.id.as_str()
            )));
        }
    }

    Ok(())
}
