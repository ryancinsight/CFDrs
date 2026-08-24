#![allow(missing_docs)]
//! Baseline CFD sanity check.
//!
//! This is the minimal end-to-end setup that exercises the three primitives
//! introduced in `docs/book/foundations.md`:
//!
//! 1. **3D incompressible flow field** — `FlowField::<f64>::new(nx, ny, nz)`
//! 2. **Diagnostic field operations** — `FlowOperations::divergence`,
//!    `FlowOperations::vorticity`, `FlowOperations::kinetic_energy`,
//!    `FlowOperations::enstrophy`
//! 3. **Reynolds-number classification** — `ReynoldsNumber::new` against
//!    a `FlowGeometry::Pipe` carrier
//!
//! Run with: `cargo run -p cfd-core --example cfd_demo`

use cfd_core::error::Result;
use cfd_core::physics::fluid_dynamics::{FlowField, FlowOperations};
use cfd_core::physics::values::{FlowGeometry, ReynoldsNumber};

fn main() -> Result<()> {
    // -------------------------------------------------------------------
    // 1. Allocate a 3D incompressible flow field. The default initial
    //    condition has zero velocity everywhere, so the divergence,
    //    vorticity, kinetic-energy, and enstrophy fields are all
    //    zero and continuity is satisfied exactly (∇·u = 0).
    // -------------------------------------------------------------------
    let nx = 32_usize;
    let ny = 32_usize;
    let nz = 32_usize;
    let flow_field = FlowField::<f64>::new(nx, ny, nz);

    let divergence = FlowOperations::divergence(&flow_field.velocity);
    let vorticity = FlowOperations::vorticity(&flow_field.velocity);
    let ke_field = FlowOperations::kinetic_energy(&flow_field.velocity);
    let enstrophy_field = FlowOperations::enstrophy(&flow_field.velocity);

    assert!(
        divergence.iter().all(|d| d.abs() < f64::EPSILON),
        "zero initial condition must satisfy continuity exactly (∇·u = 0)"
    );
    let vorticity_max: f64 = vorticity.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
    assert!(
        vorticity_max < f64::EPSILON,
        "quiescent field has zero vorticity everywhere"
    );
    let ke_max: f64 = ke_field.iter().fold(0.0_f64, |a, b| f64::max(a, *b));
    assert!(
        ke_max < f64::EPSILON,
        "quiescent field has zero kinetic energy everywhere"
    );
    let enstrophy_max: f64 = enstrophy_field.iter().fold(0.0_f64, |a, b| f64::max(a, *b));
    assert!(
        enstrophy_max < f64::EPSILON,
        "quiescent field has zero enstrophy everywhere"
    );

    // -------------------------------------------------------------------
    // 2. Reynolds-number classification. Re = 2300 sits in the
    //    transitional regime (2100 < Re < 4000) for pipe flow.
    // -------------------------------------------------------------------
    let re = ReynoldsNumber::new(2300.0, FlowGeometry::Pipe)?;
    assert!(re.is_transitional());

    Ok(())
}
