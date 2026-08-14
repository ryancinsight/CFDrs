#![expect(missing_docs, reason = "ratchet CFDRS-DOCS-1: pre-existing debt")]
#![expect(
    clippy::print_stdout,
    reason = "ratchet CFDRS-PRINT-1: pre-existing debt"
)]

use cfd_schematics::interface::presets::symmetric_bifurcation;

fn main() {
    let bp = symmetric_bifurcation("test", 0.010, 0.010, 0.004, 0.003);
    println!("Nodes:");
    for n in &bp.nodes {
        println!("  {:?}: {:?}", n.id, n.kind);
    }
    println!("Channels:");
    for c in &bp.channels {
        println!("  {:?} ({:?} -> {:?})", c.id, c.from, c.to);
    }
}
