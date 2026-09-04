//! Interface traction extracted from the flow state.
//!
//! The fluid half of the Phase 0 fluid-structure coupling (atlas ADR 0059).
//! CFDrs owns this because it reads its own pressure, velocity, and viscosity;
//! the structural solve that consumes the result lives elsewhere, and no
//! dependency runs between them.
//!
//! ## Traction on a fluid-solid interface
//!
//! For an incompressible Newtonian fluid the Cauchy stress is
//!
//! ```text
//! sigma = -p I + tau,      tau = mu (grad u + grad u^T)
//! ```
//!
//! and the traction the fluid exerts on the solid is
//!
//! ```text
//! t = sigma . n_s = -p n_s + tau . n_s
//! ```
//!
//! where `n_s` is the **solid's outward normal**, pointing from the solid into
//! the fluid. That is the Cauchy convention: `sigma . n` is the force per area
//! on a surface exerted by the material the normal points toward, and here the
//! fluid is what `n_s` points toward.
//!
//! The direction is easy to invert and the consequence is silent. Under a
//! positive pressure the correct load is compressive, pushing into the solid;
//! taking the fluid-outward normal instead yields `+p` outward, a tension that
//! no downstream structural oracle would flag as unphysical. It is fixed here
//! and asserted by [`tests::pressure_loads_the_solid_in_compression`].
//!
//! ## Interface discretisation
//!
//! The interface is the set of cell faces separating a fluid cell from a solid
//! cell, read from [`SimulationFields::mask`]. On a structured grid those faces
//! are axis-aligned, so each normal is one of `(±1, 0)` or `(0, ±1)` exactly —
//! no geometric reconstruction, and no rounding in the normal itself.
//!
//! This is the conforming-interface case atlas ADR 0059 requires: face
//! midpoints are the interface nodes, and a consumer meshing its structure on
//! the same faces has one-to-one correspondence.
//!
//! ## Units
//!
//! Values are the crate's native raw scalars, matching every other field here.
//! Pressure and traction are Pa, velocity m/s, viscosity Pa·s, spacing m.
//! Adapting these to a typed `harmonia::FieldEnvelope` is a separate step,
//! dependency-ordered after it per atlas ADR 0050, and gated on the aequitas
//! stress semantics marker so traction cannot be filled from a pressure.

use crate::fields::{Field2D, SimulationFields};
use cfd_core::CfdScalar;

/// Which face of a fluid cell the interface crosses.
///
/// A two-variant-per-axis enum rather than a signed index pair: the four
/// axis-aligned outward normals are the entire domain, and naming them keeps
/// an invalid normal unrepresentable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InterfaceFace {
    /// Solid lies at increasing x; solid outward normal `(-1, 0)`.
    East,
    /// Solid lies at decreasing x; solid outward normal `(+1, 0)`.
    West,
    /// Solid lies at increasing y; solid outward normal `(0, -1)`.
    North,
    /// Solid lies at decreasing y; solid outward normal `(0, +1)`.
    South,
}

impl InterfaceFace {
    /// The solid's outward unit normal, pointing from the solid into the fluid.
    ///
    /// This is the normal the Cauchy traction `sigma . n` is taken against, so
    /// it is the solid's, not the fluid's. `East` names where the solid sits
    /// relative to the fluid cell, so its normal points back the other way.
    ///
    /// Exact by construction on a structured grid: every component is `0`,
    /// `1`, or `-1`, so the normal contributes no rounding to the traction.
    #[must_use]
    pub fn solid_outward_normal<T: CfdScalar + Copy>(self) -> (T, T) {
        match self {
            Self::East => (-T::ONE, T::ZERO),
            Self::West => (T::ONE, T::ZERO),
            Self::North => (T::ZERO, -T::ONE),
            Self::South => (T::ZERO, T::ONE),
        }
    }

    /// Index of the solid cell across this face from fluid cell `(i, j)`.
    ///
    /// Returns `None` at a domain edge, where the neighbour does not exist.
    #[must_use]
    pub fn neighbour(self, i: usize, j: usize, nx: usize, ny: usize) -> Option<(usize, usize)> {
        match self {
            Self::East if i + 1 < nx => Some((i + 1, j)),
            Self::West if i > 0 => Some((i - 1, j)),
            Self::North if j + 1 < ny => Some((i, j + 1)),
            Self::South if j > 0 => Some((i, j - 1)),
            _ => None,
        }
    }

    /// Every face, in a fixed order.
    ///
    /// The order is part of the interface node ordering a coupling boundary
    /// marshals against, so it is stated rather than incidental.
    #[must_use]
    pub const fn all() -> [Self; 4] {
        [Self::East, Self::West, Self::North, Self::South]
    }
}

/// Traction on one interface face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceTraction<T> {
    /// Fluid cell adjacent to the face.
    pub cell: (usize, usize),
    /// Which face of that cell the interface crosses.
    pub face: InterfaceFace,
    /// Traction x-component on the solid \[Pa].
    pub tx: T,
    /// Traction y-component on the solid \[Pa].
    pub ty: T,
}

/// Extract the traction the fluid exerts on every fluid-solid interface face.
///
/// Faces are visited in cell order, and within a cell in [`InterfaceFace::all`]
/// order, so the result is a deterministic sequence a coupling boundary can
/// marshal against.
///
/// # Viscous term
///
/// `tau . n` is evaluated from the velocity gradient at the interface cell.
/// Gradients use central differences where both neighbours are fluid and a
/// one-sided difference where they are not, because differencing across the
/// interface would read velocity from inside the solid.
#[must_use]
pub fn interface_traction<T>(fields: &SimulationFields<T>, dx: T, dy: T) -> Vec<FaceTraction<T>>
where
    T: CfdScalar + Copy,
{
    let (nx, ny) = (fields.nx, fields.ny);
    let mut out = Vec::new();

    for j in 0..ny {
        for i in 0..nx {
            if !fields.mask.at(i, j) {
                continue;
            }
            for face in InterfaceFace::all() {
                let Some((si, sj)) = face.neighbour(i, j, nx, ny) else {
                    continue;
                };
                if fields.mask.at(si, sj) {
                    continue;
                }
                let (nx_c, ny_c) = face.solid_outward_normal::<T>();
                let (tx, ty) = face_traction(fields, i, j, nx_c, ny_c, dx, dy);
                out.push(FaceTraction {
                    cell: (i, j),
                    face,
                    tx,
                    ty,
                });
            }
        }
    }
    out
}

/// Traction `sigma . n_s` on the solid at one face.
fn face_traction<T>(
    fields: &SimulationFields<T>,
    i: usize,
    j: usize,
    nx_c: T,
    ny_c: T,
    dx: T,
    dy: T,
) -> (T, T)
where
    T: CfdScalar + Copy,
{
    let two = T::ONE + T::ONE;
    let p = fields.p.at(i, j);
    let mu = fields.viscosity.at(i, j);

    let du_dx = gradient_x(&fields.u, &fields.mask, i, j, dx);
    let du_dy = gradient_y(&fields.u, &fields.mask, i, j, dy);
    let dv_dx = gradient_x(&fields.v, &fields.mask, i, j, dx);
    let dv_dy = gradient_y(&fields.v, &fields.mask, i, j, dy);

    // tau = mu (grad u + grad u^T), symmetric for a Newtonian fluid.
    let tau_xx = mu * two * du_dx;
    let tau_yy = mu * two * dv_dy;
    let tau_xy = mu * (du_dy + dv_dx);

    // t = sigma . n_s, with sigma = -p I + tau
    let tx = (tau_xx - p) * nx_c + tau_xy * ny_c;
    let ty = tau_xy * nx_c + (tau_yy - p) * ny_c;
    (tx, ty)
}

/// `d/dx` of `field` at `(i, j)`, never differencing into a solid cell.
fn gradient_x<T>(field: &Field2D<T>, mask: &Field2D<bool>, i: usize, j: usize, dx: T) -> T
where
    T: CfdScalar + Copy,
{
    let two = T::ONE + T::ONE;
    let nx = field.nx();
    let east = i + 1 < nx && mask.at(i + 1, j);
    let west = i > 0 && mask.at(i - 1, j);
    match (east, west) {
        (true, true) => (field.at(i + 1, j) - field.at(i - 1, j)) / (two * dx),
        (true, false) => (field.at(i + 1, j) - field.at(i, j)) / dx,
        (false, true) => (field.at(i, j) - field.at(i - 1, j)) / dx,
        (false, false) => T::ZERO,
    }
}

/// `d/dy` of `field` at `(i, j)`, never differencing into a solid cell.
fn gradient_y<T>(field: &Field2D<T>, mask: &Field2D<bool>, i: usize, j: usize, dy: T) -> T
where
    T: CfdScalar + Copy,
{
    let two = T::ONE + T::ONE;
    let ny = field.ny();
    let north = j + 1 < ny && mask.at(i, j + 1);
    let south = j > 0 && mask.at(i, j - 1);
    match (north, south) {
        (true, true) => (field.at(i, j + 1) - field.at(i, j - 1)) / (two * dy),
        (true, false) => (field.at(i, j + 1) - field.at(i, j)) / dy,
        (false, true) => (field.at(i, j) - field.at(i, j - 1)) / dy,
        (false, false) => T::ZERO,
    }
}

#[cfg(test)]
mod tests {
    use super::{interface_traction, InterfaceFace};
    use crate::fields::SimulationFields;

    /// A 3x3 grid whose centre column at `j = 1` is solid, everything fluid.
    fn one_solid_cell(pressure: f64, u: f64, viscosity: f64) -> SimulationFields<f64> {
        let mut fields = SimulationFields::<f64>::new(3, 3);
        for j in 0..3 {
            for i in 0..3 {
                *fields.p.at_mut(i, j).expect("in bounds") = pressure;
                *fields.u.at_mut(i, j).expect("in bounds") = u;
                *fields.viscosity.at_mut(i, j).expect("in bounds") = viscosity;
                *fields.mask.at_mut(i, j).expect("in bounds") = true;
            }
        }
        *fields.mask.at_mut(1, 1).expect("in bounds") = false;
        fields
    }

    #[test]
    fn pressure_loads_the_solid_in_compression() {
        // The sign that is easy to invert and silent when wrong. A positive
        // pressure must push *into* the solid at every face, never out of it.
        let fields = one_solid_cell(100.0, 0.0, 0.0);
        let faces = interface_traction(&fields, 1.0, 1.0);
        assert_eq!(faces.len(), 4, "four fluid cells touch the solid cell");

        for face in &faces {
            let (dx, dy) = match face.face {
                InterfaceFace::East => (1.0, 0.0),
                InterfaceFace::West => (-1.0, 0.0),
                InterfaceFace::North => (0.0, 1.0),
                InterfaceFace::South => (0.0, -1.0),
            };
            // Traction must have a positive component along fluid -> solid.
            let into_solid = face.tx * dx + face.ty * dy;
            assert!(
                into_solid > 0.0,
                "{:?} at {:?} loads the solid outward ({into_solid}), not in compression",
                face.face,
                face.cell
            );
        }
    }

    #[test]
    fn hydrostatic_traction_is_exactly_minus_p_times_the_normal() {
        // Zero velocity kills the viscous term, leaving t = -p n_s exactly.
        let fields = one_solid_cell(100.0, 0.0, 0.5);
        for face in interface_traction(&fields, 1.0, 1.0) {
            let (nx, ny) = face.face.solid_outward_normal::<f64>();
            assert_eq!(face.tx, -100.0 * nx);
            assert_eq!(face.ty, -100.0 * ny);
        }
    }

    #[test]
    fn zero_pressure_and_zero_velocity_produce_no_load() {
        // The coupling must add no spurious loading of its own.
        let fields = one_solid_cell(0.0, 0.0, 1.0);
        for face in interface_traction(&fields, 1.0, 1.0) {
            assert_eq!(face.tx, 0.0);
            assert_eq!(face.ty, 0.0);
        }
    }

    #[test]
    fn uniform_flow_has_no_viscous_stress() {
        // grad u = 0 for a uniform field, so tau vanishes and only -p n_s
        // survives even at nonzero viscosity.
        let fields = one_solid_cell(10.0, 2.5, 0.9);
        for face in interface_traction(&fields, 1.0, 1.0) {
            let (nx, ny) = face.face.solid_outward_normal::<f64>();
            assert_eq!(face.tx, -10.0 * nx);
            assert_eq!(face.ty, -10.0 * ny);
        }
    }

    #[test]
    fn an_all_fluid_domain_has_no_interface() {
        let mut fields = one_solid_cell(1.0, 0.0, 0.0);
        *fields.mask.at_mut(1, 1).expect("in bounds") = true;
        assert!(interface_traction(&fields, 1.0, 1.0).is_empty());
    }

    #[test]
    fn normals_are_exact_and_opposite_across_an_axis() {
        // Exactness matters: a normal carrying rounding would put rounding
        // into every traction on the interface.
        let (ex, ey) = InterfaceFace::East.solid_outward_normal::<f64>();
        let (wx, wy) = InterfaceFace::West.solid_outward_normal::<f64>();
        assert_eq!((ex, ey), (-1.0, 0.0));
        assert_eq!((wx, wy), (1.0, 0.0));

        let (nx, ny) = InterfaceFace::North.solid_outward_normal::<f64>();
        let (sx, sy) = InterfaceFace::South.solid_outward_normal::<f64>();
        assert_eq!((nx, ny), (0.0, -1.0));
        assert_eq!((sx, sy), (0.0, 1.0));
    }

    #[test]
    fn shear_flow_produces_the_analytical_viscous_traction() {
        // u(y) = k y gives du/dy = k, so tau_xy = mu k and every other tau
        // component vanishes. On a south-facing interface (n_s = (0, +1)) the
        // traction is (tau_xy, -p).
        let (k, mu, p) = (3.0_f64, 2.0_f64, 7.0_f64);
        let mut fields = SimulationFields::<f64>::new(3, 3);
        for j in 0..3 {
            for i in 0..3 {
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "grid indices below 3 are exactly representable"
                )]
                let y = j as f64;
                *fields.u.at_mut(i, j).expect("in bounds") = k * y;
                *fields.p.at_mut(i, j).expect("in bounds") = p;
                *fields.viscosity.at_mut(i, j).expect("in bounds") = mu;
                *fields.mask.at_mut(i, j).expect("in bounds") = true;
            }
        }
        // Solid strip along the top row, so the interface at j = 1 faces north.
        for i in 0..3 {
            *fields.mask.at_mut(i, 2).expect("in bounds") = false;
        }

        let faces = interface_traction(&fields, 1.0, 1.0);
        assert_eq!(faces.len(), 3, "three fluid cells touch the solid strip");
        for face in faces {
            assert_eq!(face.face, InterfaceFace::North);
            // n_s = (0, -1): t = (tau_xy * -1, (tau_yy - p) * -1) = (-mu k, p)
            assert_eq!(face.tx, -mu * k);
            assert_eq!(face.ty, p);
        }
    }
}
