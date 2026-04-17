//! Canonical schematic preset comparison example.
//!
//! The historical name is retained for compatibility. The example now compares
//! the authoritative schematic presets that replaced the abandoned external
//! `csgrs` path.

use cfd_schematics::{serpentine_rect, venturi_rect};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let venturi = venturi_rect("compat_venturi", 1.8e-3, 0.7e-3, 0.9e-3, 2.2e-3);
    let serpentine = serpentine_rect("compat_serpentine", 4, 6.0e-3, 1.2e-3, 0.9e-3);

    venturi.validate()?;
    serpentine.validate()?;

    println!("Venturi preset geometry-authored: {}", venturi.is_geometry_authored());
    println!("Venturi channels: {}", venturi.channels.len());
    println!(
        "Venturi treatment channels: {}",
        venturi
            .topology_spec()
            .map_or(0, |spec| spec.treatment_channel_ids().len())
    );

    println!(
        "Serpentine preset geometry-authored: {}",
        serpentine.is_geometry_authored()
    );
    println!("Serpentine channels: {}", serpentine.channels.len());
    println!(
        "Serpentine has serpentine routing: {}",
        serpentine.topology_spec().is_some_and(|spec| spec.has_serpentine())
    );

    Ok(())
}
