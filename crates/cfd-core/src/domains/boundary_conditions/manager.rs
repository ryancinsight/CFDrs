//! Boundary condition manager for organizing and applying multiple conditions

use super::applicator::BoundaryConditionApplicator;
use super::geometry::BoundaryRegion;
use super::specification::BoundaryConditionSpec;
use nalgebra::RealField;
use std::collections::HashMap;

/// Boundary condition manager
pub struct BoundaryConditionManager<T: RealField + Copy> {
    /// Registered boundary regions
    regions: HashMap<String, BoundaryRegion<T>>,
    /// Registered applicators
    applicators: Vec<Box<dyn BoundaryConditionApplicator<T>>>,
}

impl<T: RealField + Copy> BoundaryConditionManager<T> {
    /// Create a new boundary condition manager
    #[must_use]
    pub fn new() -> Self {
        Self {
            regions: HashMap::new(),
            applicators: Vec::new(),
        }
    }

    /// Register a boundary region
    /// 
    /// # Errors
    /// Returns error if region with same ID already exists
    pub fn add_region(&mut self, region: BoundaryRegion<T>) -> Result<(), String> {
        if self.regions.contains_key(&region.id) {
            return Err(format!("Region '{}' already exists", region.id));
        }
        self.regions.insert(region.id.clone(), region);
        Ok(())
    }

    /// Register a boundary condition applicator
    pub fn add_applicator(&mut self, applicator: Box<dyn BoundaryConditionApplicator<T>>) {
        self.applicators.push(applicator);
    }

    /// Apply all boundary conditions to a field
    /// 
    /// # Errors
    /// Returns error if any boundary condition application fails
    pub fn apply_all(&self, field: &mut [T], time: T) -> Result<(), String> {
        for region in self.regions.values() {
            if let Some(ref condition_spec) = region.condition {
                self.apply_condition(field, condition_spec, time)?;
            }
        }
        Ok(())
    }

    /// Apply a specific boundary condition
    fn apply_condition(
        &self,
        field: &mut [T],
        spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String> {
        let condition = spec.evaluate_at_time(time);

        // Find appropriate applicator
        for applicator in &self.applicators {
            if applicator.supports(&condition) {
                return applicator.apply(field, spec, time);
            }
        }

        Err(format!(
            "No applicator found for condition type: {}",
            spec.condition_type()
        ))
    }

    /// Get a boundary region by ID
    #[must_use]
    pub fn get_region(&self, id: &str) -> Option<&BoundaryRegion<T>> {
        self.regions.get(id)
    }

    /// Update boundary condition for a region
    /// 
    /// # Errors
    /// Returns error if region with specified ID is not found
    pub fn update_condition(
        &mut self,
        region_id: &str,
        condition: BoundaryConditionSpec<T>,
    ) -> Result<(), String> {
        match self.regions.get_mut(region_id) {
            Some(region) => {
                region.condition = Some(condition);
                Ok(())
            }
            None => Err(format!("Region '{region_id}' not found")),
        }
    }

    /// Remove a boundary region
    /// 
    /// # Errors
    /// Returns error if region with specified ID is not found
    pub fn remove_region(&mut self, id: &str) -> Result<(), String> {
        if self.regions.remove(id).is_some() {
            Ok(())
        } else {
            Err(format!("Region '{id}' not found"))
        }
    }

    /// Get the number of registered regions
    #[must_use]
    pub fn num_regions(&self) -> usize {
        self.regions.len()
    }

    /// Get the number of registered applicators
    #[must_use]
    pub fn num_applicators(&self) -> usize {
        self.applicators.len()
    }

    /// Check if all regions have assigned conditions
    #[must_use]
    pub fn all_regions_have_conditions(&self) -> bool {
        self.regions
            .values()
            .all(super::geometry::BoundaryRegion::has_condition)
    }

    /// Get regions without conditions
    #[must_use]
    pub fn regions_without_conditions(&self) -> Vec<&str> {
        self.regions
            .values()
            .filter(|r| !r.has_condition())
            .map(|r| r.id.as_str())
            .collect()
    }
}

impl<T: RealField + Copy> Default for BoundaryConditionManager<T> {
    fn default() -> Self {
        Self::new()
    }
}
