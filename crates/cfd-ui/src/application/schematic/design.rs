//! Schematic design service — wraps cfd-schematics for UI interaction.

use cfd_schematics::NetworkBlueprint;

/// Create a new empty network blueprint.
#[must_use]
pub fn new_blueprint(name: &str) -> NetworkBlueprint {
    NetworkBlueprint::new(name.to_owned())
}
