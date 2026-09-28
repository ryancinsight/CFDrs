use super::*;
use crate::domain::model::{NetworkBlueprint, NodeKind};

fn count_inlets_outlets(bp: &NetworkBlueprint) -> (usize, usize) {
    let inlets = bp
        .nodes
        .iter()
        .filter(|n| n.kind == NodeKind::Inlet)
        .count();
    let outlets = bp
        .nodes
        .iter()
        .filter(|n| n.kind == NodeKind::Outlet)
        .count();
    (inlets, outlets)
}

#[test]
fn venturi_serpentine_has_single_inlet_outlet() {
    let bp = venturi_serpentine_rect("t", 2e-3, 0.5e-3, 0.5e-3, 1e-3, 6, 7.5e-3);
    let (i, o) = count_inlets_outlets(&bp);
    assert_eq!(i, 1, "must have exactly 1 inlet");
    assert_eq!(o, 1, "must have exactly 1 outlet");
}

#[test]
fn n_furcation_venturi_has_single_inlet_outlet() {
    let bp = n_furcation_venturi_rect("t", 2, 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
    let (i, o) = count_inlets_outlets(&bp);
    assert_eq!(i, 1);
    assert_eq!(o, 1);
}

#[test]
fn n_furcation_serpentine_has_single_inlet_outlet() {
    let bp = n_furcation_serpentine_rect("t", 3, 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
    let (i, o) = count_inlets_outlets(&bp);
    assert_eq!(i, 1);
    assert_eq!(o, 1);
}

#[test]
fn cell_separation_has_single_inlet_outlet() {
    let bp = cell_separation_rect("t", 22.5e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3, 22.5e-3);
    let (i, o) = count_inlets_outlets(&bp);
    assert_eq!(i, 1);
    assert_eq!(o, 1);
}
