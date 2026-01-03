//! Numerical methods service following Domain Service pattern

use super::traits::{DiscretizationScheme, LinearSystemSolver, TimeIntegrationScheme};
use super::{discretization, linear_solvers, time_integration};
use nalgebra::RealField;
use std::collections::HashMap;

/// Numerical methods service for managing and accessing numerical schemes
pub struct NumericalMethodsService<T: RealField + Copy> {
    /// Available discretization schemes
    discretization_schemes: HashMap<String, Box<dyn DiscretizationScheme<T>>>,
    /// Available time integration schemes
    time_integration_schemes: HashMap<String, Box<dyn TimeIntegrationScheme<T>>>,
    /// Available linear solvers
    linear_solvers: HashMap<String, Box<dyn LinearSystemSolver<T>>>,
}

impl<T: RealField + Copy> NumericalMethodsService<T> {
    /// Create new numerical methods service
    #[must_use]
    pub fn new() -> Self {
        let mut service = Self {
            discretization_schemes: HashMap::new(),
            time_integration_schemes: HashMap::new(),
            linear_solvers: HashMap::new(),
        };

        // Register default schemes
        service.register_discretization_scheme(
            "central".to_string(),
            Box::new(discretization::finite_difference::CentralDifference),
        );
        service.register_discretization_scheme(
            "upwind".to_string(),
            Box::new(discretization::finite_difference::UpwindDifference),
        );
        service.register_discretization_scheme(
            "downwind".to_string(),
            Box::new(discretization::finite_difference::DownwindDifference),
        );

        // Register default time integration schemes
        service.register_time_integration_scheme(
            "forward_euler".to_string(),
            Box::new(time_integration::ForwardEuler),
        );
        service.register_time_integration_scheme(
            "constant_derivative".to_string(),
            Box::new(time_integration::ConstantDerivative),
        );
        service.register_time_integration_scheme(
            "rk4".to_string(),
            Box::new(time_integration::RungeKutta4),
        );

        // Register default linear solvers
        service.register_linear_solver(
            "conjugate_gradient".to_string(),
            Box::new(linear_solvers::ConjugateGradient::default()),
        );
        service.register_linear_solver(
            "jacobi".to_string(),
            Box::new(linear_solvers::Jacobi::default()),
        );
        service.register_linear_solver(
            "gauss_seidel".to_string(),
            Box::new(linear_solvers::GaussSeidel::default()),
        );
        service
            .register_linear_solver("direct".to_string(), Box::new(linear_solvers::DirectSolver));

        service
    }

    /// Register discretization scheme
    pub fn register_discretization_scheme(
        &mut self,
        name: String,
        scheme: Box<dyn DiscretizationScheme<T>>,
    ) {
        self.discretization_schemes.insert(name, scheme);
    }

    /// Register time integration scheme
    pub fn register_time_integration_scheme(
        &mut self,
        name: String,
        scheme: Box<dyn TimeIntegrationScheme<T>>,
    ) {
        self.time_integration_schemes.insert(name, scheme);
    }

    /// Register linear solver
    pub fn register_linear_solver(&mut self, name: String, solver: Box<dyn LinearSystemSolver<T>>) {
        self.linear_solvers.insert(name, solver);
    }

    /// Get discretization scheme by name
    #[must_use]
    pub fn get_discretization_scheme(&self, name: &str) -> Option<&dyn DiscretizationScheme<T>> {
        self.discretization_schemes
            .get(name)
            .map(std::convert::AsRef::as_ref)
    }

    /// Get time integration scheme by name
    #[must_use]
    pub fn get_time_integration_scheme(&self, name: &str) -> Option<&dyn TimeIntegrationScheme<T>> {
        self.time_integration_schemes
            .get(name)
            .map(std::convert::AsRef::as_ref)
    }

    /// Get linear solver by name
    #[must_use]
    pub fn get_linear_solver(&self, name: &str) -> Option<&dyn LinearSystemSolver<T>> {
        self.linear_solvers
            .get(name)
            .map(std::convert::AsRef::as_ref)
    }

    /// List available discretization schemes
    #[must_use]
    pub fn list_discretization_schemes(&self) -> Vec<String> {
        self.discretization_schemes.keys().cloned().collect()
    }

    /// List available time integration schemes
    #[must_use]
    pub fn list_time_integration_schemes(&self) -> Vec<String> {
        self.time_integration_schemes.keys().cloned().collect()
    }

    /// List available linear solvers
    #[must_use]
    pub fn list_linear_solvers(&self) -> Vec<String> {
        self.linear_solvers.keys().cloned().collect()
    }
}

impl<T: RealField + Copy> Default for NumericalMethodsService<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_service_initialization() {
        let service = NumericalMethodsService::<f64>::new();

        // Check that default schemes are registered
        assert!(service.get_discretization_scheme("central").is_some());
        assert!(service.get_discretization_scheme("upwind").is_some());
        assert!(service
            .get_time_integration_scheme("forward_euler")
            .is_some());
        assert!(service.get_time_integration_scheme("rk4").is_some());
        assert!(service.get_linear_solver("conjugate_gradient").is_some());
        assert!(service.get_linear_solver("direct").is_some());
    }

    #[test]
    fn test_scheme_registration() {
        let mut service = NumericalMethodsService::<f64>::new();

        // Register a custom scheme
        struct CustomScheme;
        impl DiscretizationScheme<f64> for CustomScheme {
            fn discretize(&self, field: &[f64], _grid_spacing: f64) -> Vec<f64> {
                field.to_vec()
            }
            fn name(&self) -> &'static str {
                "Custom"
            }
            fn order(&self) -> usize {
                1
            }
        }

        service.register_discretization_scheme("custom".to_string(), Box::new(CustomScheme));
        assert!(service.get_discretization_scheme("custom").is_some());
    }

    #[test]
    fn test_list_schemes() {
        let service = NumericalMethodsService::<f64>::new();

        let disc_schemes = service.list_discretization_schemes();
        assert!(disc_schemes.contains(&"central".to_string()));
        assert!(disc_schemes.contains(&"upwind".to_string()));

        let time_schemes = service.list_time_integration_schemes();
        assert!(time_schemes.contains(&"forward_euler".to_string()));
        assert!(time_schemes.contains(&"rk4".to_string()));

        let solvers = service.list_linear_solvers();
        assert!(solvers.contains(&"conjugate_gradient".to_string()));
        assert!(solvers.contains(&"direct".to_string()));
    }
}
