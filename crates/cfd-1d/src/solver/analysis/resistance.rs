//! Resistance analysis module for network components.
//!
//! ## Theorem: Series/Parallel Resistance Duality
//!
//! **Theorem**: For a set of `n` positive hydraulic resistances `{R_1, …, R_n}`:
//!
//! 1. **Series**: `R_series = Σ_{i=1}^{n} R_i`, and `R_series ≥ max{R_i}`.
//! 2. **Parallel**: `R_parallel = 1 / Σ (1/R_i)`, and `R_parallel ≤ min{R_i}`.
//! 3. **Duality bound**: `R_parallel ≤ min{R_i} ≤ max{R_i} ≤ R_series`.
//!
//! **Proof of (3)**:
//! - `R_series = Σ R_i ≥ R_j` for any `j` (all terms positive). ✓
//! - `1/R_parallel = Σ (1/R_i) ≥ 1/R_j` for any `j`, so `R_parallel ≤ R_j`. ✓
//!
//! These inequalities are strict when `n > 1` and all `R_i` are finite and positive. ∎

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

/// Combine a collection of resistances in series.
///
/// The reduction is exact for any 1D chain of components arranged end-to-end.
pub(crate) fn series_resistance<T, I>(resistances: I) -> T
where
    T: RealField + Copy,
    I: IntoIterator<Item = T>,
{
    resistances
        .into_iter()
        .fold(T::zero(), |total, resistance| total + resistance)
}

/// Combine a collection of positive resistances in parallel.
///
/// Non-positive inputs are ignored to avoid division by zero and to keep the
/// reduction numerically stable on malformed inputs.
pub(crate) fn parallel_resistance<T, I>(resistances: I) -> T
where
    T: RealField + Copy,
    I: IntoIterator<Item = T>,
{
    let mut reciprocal_sum = T::zero();
    let mut valid_branch_count = 0usize;

    for resistance in resistances {
        if resistance > T::zero() {
            reciprocal_sum += T::one() / resistance;
            valid_branch_count += 1;
        }
    }

    if valid_branch_count > 0 && reciprocal_sum > T::zero() {
        T::one() / reciprocal_sum
    } else {
        T::zero()
    }
}

impl<T: RealField + Copy + Sum> ResistanceAnalysis<T> {
    /// Create a new resistance analysis
    #[must_use]
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
        self.resistances.insert(id, resistance);
        self.total_resistance += resistance;
    }

    /// Add resistance by type
    pub fn add_resistance_by_type(&mut self, component_type: String, resistance: T) {
        *self
            .resistance_by_type
            .entry(component_type)
            .or_insert(T::zero()) += resistance;
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
        parallel_resistance(self.resistances.values().copied())
    }

    /// Calculate series resistance
    pub fn series_resistance(&self) -> T {
        series_resistance(self.resistances.values().copied())
    }
}

impl<T: RealField + Copy + Sum> Default for ResistanceAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}
