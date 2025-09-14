//! Computational domain representations.

mod common;
mod domain_1d;
mod domain_2d;
mod domain_3d;
mod any_domain;

pub use common::{Domain, order};
pub use domain_1d::Domain1D;
pub use domain_2d::Domain2D;
pub use domain_3d::Domain3D;
pub use any_domain::AnyDomain;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Point1, Point2, Point3};
    use approx::assert_relative_eq;

    #[test]
    fn test_domain_1d() {
        let domain = Domain1D::new(0.0, 1.0);
        assert_eq!(domain.dimension(), 1);
        assert_relative_eq!(domain.length(), 1.0);
        assert!(domain.contains(&Point1::new(0.5)));
        assert!(!domain.contains(&Point1::new(1.5)));

        // Test automatic ordering
        let domain_reversed = Domain1D::new(1.0, 0.0);
        assert_eq!(domain_reversed.start, 0.0);
        assert_eq!(domain_reversed.end, 1.0);
        assert_relative_eq!(domain_reversed.length(), 1.0);
    }

    #[test]
    fn test_domain_2d() {
        let domain = Domain2D::from_scalars(0.0, 0.0, 2.0, 3.0);
        assert_eq!(domain.dimension(), 2);
        assert_relative_eq!(domain.area(), 6.0);
        assert!(domain.contains(&Point2::new(1.0, 1.0)));
        assert!(!domain.contains(&Point2::new(3.0, 1.0)));

        // Test from_points constructor
        let domain2 = Domain2D::from_points(Point2::new(0.0, 0.0), Point2::new(2.0, 3.0));
        assert_relative_eq!(domain2.area(), 6.0);
    }

    #[test]
    fn test_domain_3d() {
        let domain = Domain3D::new(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 3.0, 4.0));
        assert_eq!(domain.dimension(), 3);
        assert_relative_eq!(domain.volume(), 24.0);
        assert!(domain.contains(&Point3::new(1.0, 1.0, 1.0)));
        assert!(!domain.contains(&Point3::new(3.0, 1.0, 1.0)));
    }

    #[test]
    fn test_any_domain_from_conversions() {
        let domain_1d = Domain1D::new(0.0, 1.0);
        let any_domain: AnyDomain<f64> = domain_1d.into();
        assert_eq!(any_domain.dimension(), 1);

        let domain_2d = Domain2D::from_scalars(0.0, 0.0, 1.0, 1.0);
        let any_domain: AnyDomain<f64> = domain_2d.into();
        assert_eq!(any_domain.dimension(), 2);

        let domain_3d = Domain3D::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let any_domain: AnyDomain<f64> = domain_3d.into();
        assert_eq!(any_domain.dimension(), 3);
    }

    #[test]
    fn test_dimension_specific_contains_methods() {
        let domain_1d = Domain1D::new(0.0, 1.0);
        let domain_2d = Domain2D::from_scalars(0.0, 0.0, 1.0, 1.0);
        let domain_3d = Domain3D::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));

        // Test 1D domain
        assert_eq!(domain_1d.contains_1d(&Point1::new(0.5)), Some(true));
        assert_eq!(domain_1d.contains_1d(&Point1::new(1.5)), Some(false));
        assert_eq!(domain_1d.contains_2d(&Point2::new(0.5, 0.5)), None);
        assert_eq!(domain_1d.contains_3d(&Point3::new(0.5, 0.5, 0.5)), None);

        // Test 2D domain
        assert_eq!(domain_2d.contains_1d(&Point1::new(0.5)), None);
        assert_eq!(domain_2d.contains_2d(&Point2::new(0.5, 0.5)), Some(true));
        assert_eq!(domain_2d.contains_2d(&Point2::new(1.5, 0.5)), Some(false));
        assert_eq!(domain_2d.contains_3d(&Point3::new(0.5, 0.5, 0.5)), None);

        // Test 3D domain
        assert_eq!(domain_3d.contains_1d(&Point1::new(0.5)), None);
        assert_eq!(domain_3d.contains_2d(&Point2::new(0.5, 0.5)), None);
        assert_eq!(
            domain_3d.contains_3d(&Point3::new(0.5, 0.5, 0.5)),
            Some(true)
        );
        assert_eq!(
            domain_3d.contains_3d(&Point3::new(1.5, 0.5, 0.5)),
            Some(false)
        );
    }

    #[test]
    fn test_any_domain_dimension_specific_contains() {
        let any_1d: AnyDomain<f64> = Domain1D::new(0.0, 1.0).into();
        let any_2d: AnyDomain<f64> = Domain2D::from_scalars(0.0, 0.0, 1.0, 1.0).into();
        let any_3d: AnyDomain<f64> =
            Domain3D::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)).into();

        // Test proper dimension matches
        assert_eq!(any_1d.contains_1d(&Point1::new(0.5)), Some(true));
        assert_eq!(any_2d.contains_2d(&Point2::new(0.5, 0.5)), Some(true));
        assert_eq!(any_3d.contains_3d(&Point3::new(0.5, 0.5, 0.5)), Some(true));

        // Test dimension mismatches return None
        assert_eq!(any_1d.contains_2d(&Point2::new(0.5, 0.5)), None);
        assert_eq!(any_1d.contains_3d(&Point3::new(0.5, 0.5, 0.5)), None);
        assert_eq!(any_2d.contains_1d(&Point1::new(0.5)), None);
        assert_eq!(any_2d.contains_3d(&Point3::new(0.5, 0.5, 0.5)), None);
        assert_eq!(any_3d.contains_1d(&Point1::new(0.5)), None);
        assert_eq!(any_3d.contains_2d(&Point2::new(0.5, 0.5)), None);
    }
}