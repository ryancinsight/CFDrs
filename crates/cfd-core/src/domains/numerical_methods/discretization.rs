//! Discretization schemes for spatial derivatives

use super::traits::DiscretizationScheme;
use nalgebra::RealField;

/// Finite difference schemes
pub mod finite_difference {
    use super::{DiscretizationScheme, RealField};

    /// Central difference scheme (2nd order)
    #[derive(Debug, Clone)]
    pub struct CentralDifference;

    impl<T: RealField + Copy> DiscretizationScheme<T> for CentralDifference {
        fn discretize(&self, field: &[T], grid_spacing: T) -> Vec<T> {
            if field.len() < 3 {
                return field.to_vec();
            }

            let two = T::one() + T::one();
            let dx = grid_spacing;

            field
                .windows(3)
                .map(|window| (window[2] - window[0]) / (two * dx))
                .collect()
        }

        fn name(&self) -> &'static str {
            "Central Difference"
        }

        fn order(&self) -> usize {
            2
        }
    }

    /// Upwind difference scheme (1st order)
    #[derive(Debug, Clone)]
    pub struct UpwindDifference;

    impl<T: RealField + Copy> DiscretizationScheme<T> for UpwindDifference {
        fn discretize(&self, field: &[T], grid_spacing: T) -> Vec<T> {
            if field.len() < 2 {
                return field.to_vec();
            }

            let dx = grid_spacing;

            field
                .windows(2)
                .map(|window| (window[1] - window[0]) / dx)
                .collect()
        }

        fn name(&self) -> &'static str {
            "Upwind Difference"
        }

        fn order(&self) -> usize {
            1
        }
    }

    /// Downwind difference scheme (1st order)
    #[derive(Debug, Clone)]
    pub struct DownwindDifference;

    impl<T: RealField + Copy> DiscretizationScheme<T> for DownwindDifference {
        fn discretize(&self, field: &[T], grid_spacing: T) -> Vec<T> {
            if field.len() < 2 {
                return field.to_vec();
            }

            let dx = grid_spacing;

            field
                .windows(2)
                .skip(1)
                .map(|window| (window[1] - window[0]) / dx)
                .collect()
        }

        fn name(&self) -> &'static str {
            "Downwind Difference"
        }

        fn order(&self) -> usize {
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::finite_difference::{CentralDifference, DownwindDifference, UpwindDifference};
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_central_difference_scheme() {
        let scheme = finite_difference::CentralDifference;
        let field = vec![1.0f64, 2.0, 4.0, 7.0, 11.0];
        let dx = 1.0f64;

        let result = scheme.discretize(&field, dx);

        // Central difference: (f[i+1] - f[i-1])/(2*dx)
        // Expected: [(4-1)/2, (7-2)/2, (11-4)/2] = [1.5, 2.5, 3.5]
        assert_eq!(result.len(), 3);
        assert_relative_eq!(result[0], 1.5, epsilon = 1e-10);
        assert_relative_eq!(result[1], 2.5, epsilon = 1e-10);
        assert_relative_eq!(result[2], 3.5, epsilon = 1e-10);

        assert_eq!(
            <CentralDifference as DiscretizationScheme<f64>>::name(&scheme),
            "Central Difference"
        );
        assert_eq!(
            <CentralDifference as DiscretizationScheme<f64>>::order(&scheme),
            2
        );
    }

    #[test]
    fn test_upwind_difference_scheme() {
        let scheme = finite_difference::UpwindDifference;
        let field = vec![1.0f64, 2.0, 4.0, 7.0];
        let dx = 1.0f64;

        let result = scheme.discretize(&field, dx);

        // Upwind difference: (f[i+1] - f[i])/dx
        // Expected: [(2-1)/1, (4-2)/1, (7-4)/1] = [1.0, 2.0, 3.0]
        assert_eq!(result.len(), 3);
        assert_relative_eq!(result[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(result[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(result[2], 3.0, epsilon = 1e-10);

        assert_eq!(
            <UpwindDifference as DiscretizationScheme<f64>>::name(&scheme),
            "Upwind Difference"
        );
        assert_eq!(
            <UpwindDifference as DiscretizationScheme<f64>>::order(&scheme),
            1
        );
    }

    #[test]
    fn test_downwind_difference_scheme() {
        let scheme = finite_difference::DownwindDifference;
        let field = vec![1.0f64, 3.0, 6.0, 10.0];
        let dx = 1.0f64;

        let result = scheme.discretize(&field, dx);

        // Downwind difference uses forward values
        assert_eq!(result.len(), 2);
        assert_relative_eq!(result[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(result[1], 4.0, epsilon = 1e-10);

        assert_eq!(
            <DownwindDifference as DiscretizationScheme<f64>>::name(&scheme),
            "Downwind Difference"
        );
        assert_eq!(
            <DownwindDifference as DiscretizationScheme<f64>>::order(&scheme),
            1
        );
    }
}
