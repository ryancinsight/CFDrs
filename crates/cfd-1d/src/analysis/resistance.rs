//! Resistance analysis module for network components.

use nalgebra::RealField;
use std::collections::HashMap;
use std::iter::Sum;

/// Resistance analysis for network components
#[derive(Debug, Clone)]
pub struct ResistanceAnalysis<T: RealField + Copy> {
    /// Hydraulic resistances [Pa·s/m³]
    pub resistances: HashMap<String, T>,
    /// Equivalent circuit resistance [Pa·s/m³]
    pub total_resistance: T,
    /// Resistance contributions by component type
    pub resistance_by_type: HashMap<String, T>,
    /// Critical resistance paths
    pub critical_paths: Vec<Vec<String>>,
}

impl<T: RealField + Copy + Sum> ResistanceAnalysis<T> {
    /// Create a new resistance analysis
    pub fn new() -> Self {
        Self {
            resistances: HashMap::new(),
            total_resistance: T::zero(),
            resistance_by_type: HashMap::new(),
            critical_paths: Vec::new(),
        }
    }

    /// Add resistance data for a component
    pub fn add_resistance(&mut self, id: String, resistance: T) {
        self.resistances.insert(id.clone(), resistance);
        self.total_resistance = self.total_resistance + resistance;
    }

    /// Add resistance by type
    pub fn add_resistance_by_type(&mut self, component_type: String, resistance: T) {
        *self.resistance_by_type.entry(component_type).or_insert(T::zero()) += resistance;
    }

    /// Add a critical path
    pub fn add_critical_path(&mut self, path: Vec<String>) {
        self.critical_paths.push(path);
    }

    /// Get the average resistance
    pub fn average_resistance(&self) -> T {
        if self.resistances.is_empty() {
            T::zero()
        } else {
            let sum: T = self.resistances.values().copied().sum();
            sum / T::from_usize(self.resistances.len()).unwrap_or_else(T::one)
        }
    }

    /// Get the maximum resistance component
    pub fn max_resistance(&self) -> Option<(&String, T)> {
        self.resistances
            .iter()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(id, &r)| (id, r))
    }

    /// Get the minimum resistance component
    pub fn min_resistance(&self) -> Option<(&String, T)> {
        self.resistances
            .iter()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(id, &r)| (id, r))
    }

    /// Calculate parallel resistance
    pub fn parallel_resistance(&self) -> T {
        if self.resistances.is_empty() {
            T::zero()
        } else {
            let sum_inv: T = self.resistances
                .values()
                .filter(|&&r| r > T::zero())
                .map(|&r| T::one() / r)
                .sum();
            
            if sum_inv > T::zero() {
                T::one() / sum_inv
            } else {
                T::zero()
            }
        }
    }

    /// Calculate series resistance
    pub fn series_resistance(&self) -> T {
        self.resistances.values().copied().sum()
    }
}

impl<T: RealField + Copy + Sum> Default for ResistanceAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}