//! Dumps the nodes and channels of a symmetric bifurcation preset.

use std::io::{self, Write};

use cfd_schematics::interface::presets::symmetric_bifurcation;

fn main() -> io::Result<()> {
    let bp = symmetric_bifurcation("test", 0.010, 0.010, 0.004, 0.003);
    let standard_output = io::stdout();
    let mut terminal = standard_output.lock();
    writeln!(terminal, "Nodes:")?;
    for n in &bp.nodes {
        writeln!(terminal, "  {:?}: {:?}", n.id, n.kind)?;
    }
    writeln!(terminal, "Channels:")?;
    for c in &bp.channels {
        writeln!(terminal, "  {:?} ({:?} -> {:?})", c.id, c.from, c.to)?;
    }
    Ok(())
}
