//! Boundary conditions domain - Physical constraints and boundary condition management.
//!
//! This module encapsulates boundary condition knowledge following DDD principles.
//! It provides abstractions for different boundary condition types and their application.

use crate::boundary::BoundaryCondition;
use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Boundary condition specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryConditionSpec<T: RealField + Copy> {
    /// Boundary condition
    pub condition: BoundaryCondition<T>,
    /// Boundary region identifier
    pub region_id: String,
    /// Time-dependent specification
    pub time_dependent: Option<TimeDependentSpec<T>>,
}

/// Time-dependent boundary condition specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeDependentSpec<T: RealField + Copy> {
    /// Time function type
    pub function_type: TimeFunctionType,
    /// Function parameters
    pub parameters: Vec<T>,
}

/// Time function types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum TimeFunctionType {
    /// Constant value
    Constant,
    /// Linear ramp
    Linear,
    /// Sinusoidal variation
    Sinusoidal,
    /// Exponential decay/growth
    Exponential,
    /// Custom function (user-defined)
    Custom(String),
}

/// Boundary condition applicator abstraction
pub trait BoundaryConditionApplicator<T: RealField + Copy>: Send + Sync {
    /// Apply boundary condition to field
    fn apply(
        &self,
        field: &mut [T],
        boundary_spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String>;
    
    /// Get applicator name
    fn name(&self) -> &str;
    
    /// Check if this applicator supports the given boundary condition
    fn supports(&self, condition: &BoundaryCondition<T>) -> bool;
}

/// Boundary region specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryRegion<T: RealField + Copy> {
    /// Region identifier
    pub id: String,
    /// Region geometry
    pub geometry: BoundaryGeometry<T>,
    /// Associated boundary condition
    pub condition: Option<BoundaryConditionSpec<T>>,
}

/// Boundary geometry types
///
/// Defines the geometric representation of boundary regions for different dimensionalities.
/// Used to specify where boundary conditions should be applied in the computational domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BoundaryGeometry<T: RealField + Copy> {
    /// Point boundary (0D) - single point in space
    Point(Point3<T>),
    /// Line boundary (1D) - line segment defined by start and end points
    Line {
        /// Starting point of the line segment
        start: Point3<T>,
        /// Ending point of the line segment
        end: Point3<T>,
    },
    /// Surface boundary (2D) - polygonal surface defined by vertices
    Surface {
        /// Vertices defining the boundary surface (ordered)
        vertices: Vec<Point3<T>>,
    },
    /// Volume boundary (3D) - volumetric region defined by vertices
    Volume {
        /// Vertices defining the boundary volume
        vertices: Vec<Point3<T>>,
    },
}

/// Standard boundary condition applicators
pub mod applicators {
    use super::{BoundaryCondition, BoundaryConditionApplicator, BoundaryConditionSpec};
    use nalgebra::RealField;
    
    /// Dirichlet boundary condition applicator
    #[derive(Debug, Clone)]
    pub struct DirichletApplicator;
    
    impl<T: RealField + Copy> BoundaryConditionApplicator<T> for DirichletApplicator {
        fn apply(
            &self,
            field: &mut [T],
            boundary_spec: &BoundaryConditionSpec<T>,
            _time: T,
        ) -> Result<(), String> {
            match &boundary_spec.condition {
                BoundaryCondition::Dirichlet { value } => {
                    // Apply scalar Dirichlet condition
                    field.iter_mut().for_each(|f| *f = *value);
                    Ok(())
                }
                _ => Err(format!("DirichletApplicator does not support {:?}", boundary_spec.condition)),
            }
        }
        
        fn name(&self) -> &str {
            "DirichletApplicator"
        }
        
        fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
            matches!(condition, BoundaryCondition::Dirichlet { .. })
        }
    }
    
    /// Neumann boundary condition applicator
    #[derive(Debug, Clone)]
    pub struct NeumannApplicator;
    
    impl<T: RealField + Copy> BoundaryConditionApplicator<T> for NeumannApplicator {
        fn apply(
            &self,
            _field: &mut [T],
            boundary_spec: &BoundaryConditionSpec<T>,
            _time: T,
        ) -> Result<(), String> {
            match &boundary_spec.condition {
                BoundaryCondition::Neumann { flux: _ } => {
                    // Neumann conditions typically modify gradient calculations
                    // Implementation depends on discretization scheme
                    Ok(())
                }
                _ => Err(format!("NeumannApplicator does not support {:?}", boundary_spec.condition)),
            }
        }
        
        fn name(&self) -> &str {
            "NeumannApplicator"
        }
        
        fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
            matches!(condition, BoundaryCondition::Neumann { .. })
        }
    }
}

/// Time-dependent evaluator for boundary conditions
pub struct TimeDependentEvaluator<T: RealField + Copy> {
    /// Registered time functions
    functions: HashMap<String, Box<dyn Fn(T) -> T + Send + Sync>>,
}

impl<T: RealField + Copy + num_traits::Float> TimeDependentEvaluator<T> {
    /// Create a new time-dependent evaluator
    pub fn new() -> Self {
        let mut evaluator = Self {
            functions: HashMap::new(),
        };
        evaluator.register_common_functions();
        evaluator
    }
    
    /// Register a time function
    pub fn register_function<F>(&mut self, name: String, function: F)
    where
        F: Fn(T) -> T + Send + Sync + 'static,
    {
        self.functions.insert(name, Box::new(function));
    }
    
    /// Evaluate a registered function
    pub fn evaluate(&self, name: &str, time: T) -> Option<T> {
        self.functions.get(name).map(|f| f(time))
    }
    
    /// Register common time-dependent functions
    pub fn register_common_functions(&mut self) {
        // Sinusoidal function: sin(ωt)
        self.register_function("sin".to_string(), |t| num_traits::Float::sin(t));
        
        // Exponential decay: exp(-t)
        self.register_function("exp_decay".to_string(), |t| num_traits::Float::exp(-t));
        
        // Step function: H(t-1) where H is Heaviside
        self.register_function("step".to_string(), |t| {
            if t >= T::one() {
                T::one()
            } else {
                T::zero()
            }
        });
        
        // Ramp function: max(0, t)
        self.register_function("ramp".to_string(), |t| {
            if t >= T::zero() {
                t
            } else {
                T::zero()
            }
        });
    }
    
    /// Evaluate time-dependent specification
    pub fn evaluate_spec(&mut self, spec: &TimeDependentSpec<T>, time: T) -> T {
        match spec.function_type {
            TimeFunctionType::Constant => {
                if !spec.parameters.is_empty() {
                    spec.parameters[0]
                } else {
                    T::one()
                }
            }
            TimeFunctionType::Linear => {
                // Linear: a*t + b
                let a = spec.parameters.get(0).copied().unwrap_or(T::one());
                let b = spec.parameters.get(1).copied().unwrap_or(T::zero());
                a * time + b
            }
            TimeFunctionType::Sinusoidal => {
                // Sinusoidal: A*sin(ωt + φ) + offset
                let amplitude = spec.parameters.get(0).copied().unwrap_or(T::one());
                let frequency = spec.parameters.get(1).copied().unwrap_or(T::one());
                let phase = spec.parameters.get(2).copied().unwrap_or(T::zero());
                let offset = spec.parameters.get(3).copied().unwrap_or(T::zero());
                amplitude * num_traits::Float::sin(frequency * time + phase) + offset
            }
            TimeFunctionType::Exponential => {
                // Exponential: A*exp(λt) + offset
                let amplitude = spec.parameters.get(0).copied().unwrap_or(T::one());
                let rate = spec.parameters.get(1).copied().unwrap_or(-T::one());
                let offset = spec.parameters.get(2).copied().unwrap_or(T::zero());
                amplitude * num_traits::Float::exp(rate * time) + offset
            }
            TimeFunctionType::Custom(ref name) => {
                self.evaluate(name, time).unwrap_or(T::one())
            }
        }
    }
}

impl<T: RealField + Copy + num_traits::Float> Default for TimeDependentEvaluator<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Boundary conditions service for managing multiple boundary regions
pub struct BoundaryConditionsService<T: RealField + Copy> {
    /// Registered boundary regions
    regions: HashMap<String, BoundaryRegion<T>>,
    /// Registered applicators
    applicators: Vec<Box<dyn BoundaryConditionApplicator<T>>>,
}

impl<T: RealField + Copy> BoundaryConditionsService<T> {
    /// Create a new boundary conditions service
    pub fn new() -> Self {
        use applicators::{DirichletApplicator, NeumannApplicator};
        
        Self {
            regions: HashMap::new(),
            applicators: vec![
                Box::new(DirichletApplicator),
                Box::new(NeumannApplicator),
            ],
        }
    }
    
    /// Register a boundary region
    pub fn register_region(&mut self, region: BoundaryRegion<T>) {
        self.regions.insert(region.id.clone(), region);
    }
    
    /// Apply boundary conditions to field
    pub fn apply_conditions(&self, field: &mut [T], time: T) -> Result<(), String> {
        for region in self.regions.values() {
            if let Some(spec) = &region.condition {
                let applicator = self
                    .applicators
                    .iter()
                    .find(|app| app.supports(&spec.condition));
                
                if let Some(app) = applicator {
                    app.apply(field, spec, time)?;
                }
            }
        }
        Ok(())
    }
    
    /// Get boundary region by ID
    pub fn get_region(&self, id: &str) -> Option<&BoundaryRegion<T>> {
        self.regions.get(id)
    }
    
    /// List all boundary regions
    pub fn list_regions(&self) -> Vec<&str> {
        self.regions
            .keys()
            .map(std::string::String::as_str)
            .collect()
    }
}

impl<T: RealField + Copy> Default for BoundaryConditionsService<T> {
    fn default() -> Self {
        Self::new()
    }
}
