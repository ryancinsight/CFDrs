//! Boundary conditions domain - Physical constraints and boundary condition management.
//!
//! This module encapsulates boundary condition knowledge following DDD principles.
//! It provides abstractions for different boundary condition types and their application.

use crate::boundary::BoundaryCondition;
use cfd_core::numeric;
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
pub struct TimeDependentSpec<T: RealField + Copy> {
    /// Time function type
    pub function_type: TimeFunctionType,
    /// Function parameters
    pub parameters: Vec<T>,
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
/// Boundary region specification
pub struct BoundaryRegion<T: RealField + Copy> {
    /// Region identifier
    pub id: String,
    /// Region geometry
    pub geometry: BoundaryGeometry<T>,
    /// Associated boundary condition
    pub condition: Option<BoundaryConditionSpec<T>>,
/// Boundary geometry types
///
/// Defines the geometric representation of boundary regions for different dimensionalities.
/// Used to specify where boundary conditions should be applied in the computational domain.
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
    /// Volume boundary (3D) - volumetric region defined by vertices
    Volume {
        /// Vertices defining the boundary volume
/// Standard boundary condition applicators
pub mod applicators {
    use super::{BoundaryCondition, BoundaryConditionApplicator, BoundaryConditionSpec, RealField};
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
                _ => Err(
                    "Dirichlet applicator only supports Dirichlet boundary conditions".to_string(),
                ),
            }
        }
        fn name(&self) -> &str {
            "Dirichlet"
        fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
            matches!(condition, BoundaryCondition::Dirichlet { .. })
    }
    /// Neumann boundary condition applicator
    pub struct NeumannApplicator;
    impl<T: RealField + Copy> BoundaryConditionApplicator<T> for NeumannApplicator {
                BoundaryCondition::Neumann { gradient } => {
                    // Apply scalar Neumann condition
                    if let Some(last) = field.last_mut() {
                        *last = *gradient;
                    }
                _ => {
                    Err("Neumann applicator only supports Neumann boundary conditions".to_string())
            "Neumann"
            matches!(condition, BoundaryCondition::Neumann { .. })
    /// Wall boundary condition applicator
    pub struct WallApplicator;
    impl<T: RealField + Copy> BoundaryConditionApplicator<T> for WallApplicator {
            _boundary_spec: &BoundaryConditionSpec<T>,
            // Apply no-slip wall condition (zero velocity)
            field.iter_mut().for_each(|f| *f = T::zero());
            Ok(())
            "Wall"
            matches!(condition, BoundaryCondition::Wall { .. })
/// Time-dependent boundary condition evaluator
/// Provides functionality for evaluating time-varying boundary conditions
pub struct TimeDependentEvaluator<T: RealField + Copy + num_traits::Float> {
    /// Function registry for time-dependent evaluations
    functions: HashMap<String, Box<dyn Fn(T) -> T + Send + Sync>>,
impl<T: RealField + Copy + num_traits::Float> TimeDependentEvaluator<T> {
    /// Create new evaluator
    #[must_use]
    pub fn new() -> Self {
        Self {
            functions: HashMap::new(),
    /// Register a time-dependent function
    pub fn register_function<F>(&mut self, name: String, func: F)
    where
        F: Fn(T) -> T + Send + Sync + 'static,
    {
        self.functions.insert(name, Box::new(func));
    /// Evaluate a registered function at given time
    pub fn evaluate(&self, name: &str, time: T) -> Option<T> {
        self.functions.get(name).map(|f| f(time))
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
        });
        // Ramp function: max(0, t)
        self.register_function("ramp".to_string(), |t| {
            if t >= T::zero() {
                t
    /// Evaluate time-dependent specification
    pub fn evaluate_spec(&mut self, spec: &TimeDependentSpec<T>, time: T) -> T {
        match spec.function_type {
            TimeFunctionType::Constant => spec.parameters.first().copied().unwrap_or_else(T::zero),
            TimeFunctionType::Linear => {
                if spec.parameters.len() >= 2 {
                    spec.parameters[0] + spec.parameters[1] * time
                } else {
                    T::zero()
            TimeFunctionType::Sinusoidal => {
                if spec.parameters.len() >= 3 {
                    // A * sin(ω * t + φ)
                    let amplitude = spec.parameters[0];
                    let omega = spec.parameters[1];
                    let phase = spec.parameters[2];
                    // Use proper sine function
                    let angle = omega * time + phase;
                    // Convert to f64 for trig calculation, then back
                    let angle_f64: f64 = angle.to_subset().unwrap_or(0.0);
                    let sin_value = angle_f64.sin();
                    amplitude * cfd_core::numeric::from_f64(sin_value)?
            TimeFunctionType::Exponential => {
                    // A * exp(λ * t)
                    let lambda = spec.parameters[1];
                    // Use proper exponential function
                    let exponent = lambda * time;
                    // Convert to f64 for exp calculation, then back
                    let exp_f64: f64 = exponent.to_subset().unwrap_or(0.0);
                    let exp_value = exp_f64.exp();
                    amplitude * cfd_core::numeric::from_f64(exp_value)?
            TimeFunctionType::Custom(_) => {
                // Custom functions would be implemented via plugin system
/// Boundary conditions service following Domain Service pattern
pub struct BoundaryConditionsService<T: RealField + Copy> {
    /// Available applicators
    applicators: HashMap<String, Box<dyn BoundaryConditionApplicator<T>>>,
    /// Boundary regions
    regions: HashMap<String, BoundaryRegion<T>>,
impl<T: RealField + Copy> BoundaryConditionsService<T> {
    /// Create new boundary conditions service
        let mut service = Self {
            applicators: HashMap::new(),
            regions: HashMap::new(),
        };
        // Register default applicators
        service.register_applicator(
            "dirichlet".to_string(),
            Box::new(applicators::DirichletApplicator),
        );
            "neumann".to_string(),
            Box::new(applicators::NeumannApplicator),
        service.register_applicator("wall".to_string(), Box::new(applicators::WallApplicator));
        service
    /// Register boundary condition applicator
    pub fn register_applicator(
        &mut self,
        name: String,
        applicator: Box<dyn BoundaryConditionApplicator<T>>,
    ) {
        self.applicators.insert(name, applicator);
    /// Add boundary region
    pub fn add_region(&mut self, region: BoundaryRegion<T>) {
        self.regions.insert(region.id.clone(), region);
    /// Apply boundary conditions to field
    pub fn apply_boundary_conditions(&mut self, field: &mut [T], time: T) -> Result<(), String> {
        for region in self.regions.values() {
            if let Some(spec) = &region.condition {
                // Find the appropriate applicator based on the boundary condition type
                let applicator = self
                    .applicators
                    .values()
                    .find(|app| app.supports(&spec.condition));
                if let Some(app) = applicator {
                    app.apply(field, spec, time)?;
        Ok(())
    /// Get boundary region by ID
    pub fn get_region(&self, id: &str) -> Option<&BoundaryRegion<T>> {
        self.regions.get(id)
    /// List all boundary regions
    pub fn list_regions(&self) -> Vec<&str> {
        self.regions
            .keys()
            .map(std::string::String::as_str)
            .collect()
impl<T: RealField + Copy> Default for BoundaryConditionsService<T> {
    fn default() -> Self {
        Self::new()
impl<T: RealField + Copy + num_traits::Float> Default for TimeDependentEvaluator<T> {
