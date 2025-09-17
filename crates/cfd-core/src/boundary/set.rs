//! Boundary condition collection management

use super::BoundaryCondition;
use crate::error::{Error, Result};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Collection of boundary conditions
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundaryConditionSet<T: RealField + Copy> {
    conditions: HashMap<String, BoundaryCondition<T>>,
}

impl<T: RealField + Copy> BoundaryConditionSet<T> {
    /// Create empty set
    #[must_use]
    pub fn new() -> Self {
        Self {
            conditions: HashMap::new(),
        }
    }

    /// Add boundary condition
    pub fn add(&mut self, name: impl Into<String>, condition: BoundaryCondition<T>) -> &mut Self {
        self.conditions.insert(name.into(), condition);
        self
    }

    /// Get boundary condition by name
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&BoundaryCondition<T>> {
        self.conditions.get(name)
    }

    /// Get number of boundary conditions
    #[must_use]
    pub fn len(&self) -> usize {
        self.conditions.len()
    }

    /// Check if set is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.conditions.is_empty()
    }

    /// Get mutable boundary condition
    pub fn get_mut(&mut self, name: &str) -> Option<&mut BoundaryCondition<T>> {
        self.conditions.get_mut(name)
    }

    /// Remove boundary condition
    pub fn remove(&mut self, name: &str) -> Option<BoundaryCondition<T>> {
        self.conditions.remove(name)
    }

    /// Check if set contains boundary
    pub fn contains(&self, name: &str) -> bool {
        self.conditions.contains_key(name)
    }

    /// Iterate over boundaries
    pub fn iter(&self) -> impl Iterator<Item = (&String, &BoundaryCondition<T>)> {
        self.conditions.iter()
    }

    /// Validate periodic boundary pairs
    ///
    /// # Errors
    /// Returns error if periodic boundary references non-existent partner
    /// or if partner boundary is not also periodic
    pub fn validate_periodic(&self) -> Result<()> {
        for (name, bc) in &self.conditions {
            if let BoundaryCondition::Periodic { partner } = bc {
                // Check partner exists
                let partner_bc = self.conditions.get(partner).ok_or_else(|| {
                    Error::InvalidConfiguration(format!(
                        "Periodic boundary '{name}' references non-existent partner '{partner}'"
                    ))
                })?;

                // Check partner is also periodic
                if let BoundaryCondition::Periodic {
                    partner: partner_ref,
                } = partner_bc
                {
                    // Check mutual reference
                    if partner_ref != name {
                        return Err(Error::InvalidConfiguration(format!(
                            "Periodic boundaries '{name}' and '{partner}' do not reference each other"
                        )));
                    }
                } else {
                    return Err(Error::InvalidConfiguration(format!(
                        "Periodic boundary '{name}' partner '{partner}' is not periodic"
                    )));
                }
            }
        }
        Ok(())
    }
}

impl<T: RealField + Copy> Default for BoundaryConditionSet<T> {
    fn default() -> Self {
        Self::new()
    }
}
