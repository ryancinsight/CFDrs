//! `geometry/generator/shell` — Shell cuboid constructor.
//!
//! Provides a single free function [`create_shell_cuboid`] that validates
//! parameters and delegates to [`ShellCuboid::new`].

use aequitas::systems::si::quantities::Length;

use crate::error::GeometryResult;
use crate::geometry::types::ShellCuboid;

/// Create a shell millifluidic cuboid schematic.
///
/// The resulting [`ShellCuboid`] has no internal channels — the interior is
/// a single open cavity bounded by a uniform wall of `shell_thickness` on
/// all four sides. Inlet and outlet ports are represented as short line stubs
/// on the left and right walls at the horizontal midline.
///
/// # Arguments
///
/// * `outer_dims` — outer bounding-box lengths stored in base metres.
/// * `shell_thickness` — uniform wall thickness stored in base metres.
///
/// # Errors
///
/// Returns [`GeometryError`](crate::error::GeometryError) when any input
/// violates the geometric invariants (see [`ShellCuboid::new`]).
///
/// # Examples
///
/// ```rust
/// use aequitas::systems::si::quantities::Length;
/// use aequitas::systems::si::units::Millimeter;
/// use cfd_schematics::geometry::generator::shell::create_shell_cuboid;
///
/// let sc = create_shell_cuboid(
///     (
///         Length::from_unit::<Millimeter>(80.0),
///         Length::from_unit::<Millimeter>(40.0),
///     ),
///     Length::from_unit::<Millimeter>(2.0),
/// )
/// .expect("valid shell cuboid");
/// assert_eq!(sc.box_outline.len(), 10);
/// ```
pub fn create_shell_cuboid(
    outer_dims: (Length<f64>, Length<f64>),
    shell_thickness: Length<f64>,
) -> GeometryResult<ShellCuboid> {
    ShellCuboid::new(outer_dims, shell_thickness)
}

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::units::Millimeter;

    #[test]
    fn create_shell_cuboid_standard_case() {
        let sc = create_shell_cuboid(
            (
                Length::from_unit::<Millimeter>(80.0),
                Length::from_unit::<Millimeter>(40.0),
            ),
            Length::from_unit::<Millimeter>(2.0),
        )
        .expect("valid shell cuboid");
        assert_eq!(sc.outer_dims_mm(), (80.0, 40.0));
        assert!((sc.inner_dims_mm().0 - 76.0).abs() < 1e-9);
        assert!((sc.inner_dims_mm().1 - 36.0).abs() < 1e-9);
        assert_eq!(sc.box_outline.len(), 10);
    }

    #[test]
    fn create_shell_cuboid_rejects_impossible_thickness() {
        assert!(create_shell_cuboid(
            (
                Length::from_unit::<Millimeter>(10.0),
                Length::from_unit::<Millimeter>(4.0),
            ),
            Length::from_unit::<Millimeter>(3.0),
        )
        .is_err());
    }

    #[test]
    fn create_shell_cuboid_interchange_has_two_ports() {
        let sc = create_shell_cuboid(
            (
                Length::from_unit::<Millimeter>(80.0),
                Length::from_unit::<Millimeter>(40.0),
            ),
            Length::from_unit::<Millimeter>(2.0),
        )
        .expect("structural invariant");
        let ix = sc.to_interchange();
        assert_eq!(ix.ports.len(), 2);
        assert_eq!(ix.ports[0].label, "inlet");
        assert_eq!(ix.ports[1].label, "outlet");
        assert!((ix.shell_thickness_mm - 2.0).abs() < 1e-9);
    }
}
