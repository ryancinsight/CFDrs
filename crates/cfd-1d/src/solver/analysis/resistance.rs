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

use aequitas::systems::si::quantities::HydraulicResistance;
use cfd_core::conversion::SafeFromUsize;
use cfd_core::CfdScalar;
use std::collections::HashMap;
use std::iter::Sum;

/// Resistance analysis for network components
#[derive(Debug, Clone)]
pub struct ResistanceAnalysis<T: CfdScalar + Copy> {
    /// Hydraulic resistances [Pa·s/m³]
    pub resistances: HashMap<String, HydraulicResistance<T>>,
    /// Equivalent circuit resistance [Pa·s/m³]
    pub total_resistance: HydraulicResistance<T>,
    /// Resistance contributions by component type
    pub resistance_by_type: HashMap<String, HydraulicResistance<T>>,
    /// Critical resistance paths
    pub critical_paths: Vec<Vec<String>>,
}

/// Combine a collection of resistances in series.
///
/// The reduction is exact for any 1D chain of components arranged end-to-end.
pub(crate) fn series_resistance<T, I>(resistances: I) -> HydraulicResistance<T>
where
    T: CfdScalar + Copy,
    I: IntoIterator<Item = HydraulicResistance<T>>,
{
    HydraulicResistance::from_base(
        resistances
            .into_iter()
            .map(|resistance| resistance.into_base())
            .fold(T::ZERO, |total, resistance| total + resistance),
    )
}

/// Combine a collection of positive resistances in parallel.
///
/// Non-positive inputs are ignored to avoid division by zero and to keep the
/// reduction numerically stable on malformed inputs.
pub(crate) fn parallel_resistance<T, I>(resistances: I) -> HydraulicResistance<T>
where
    T: CfdScalar + Copy,
    I: IntoIterator<Item = HydraulicResistance<T>>,
{
    let mut reciprocal_sum = T::ZERO;
    let mut valid_branch_count = 0usize;

    for resistance in resistances {
        let resistance = resistance.into_base();
        if resistance > T::ZERO {
            reciprocal_sum += T::ONE / resistance;
            valid_branch_count += 1;
        }
    }

    if valid_branch_count > 0 && reciprocal_sum > T::ZERO {
        HydraulicResistance::from_base(T::ONE / reciprocal_sum)
    } else {
        HydraulicResistance::from_base(T::ZERO)
    }
}

impl<T: CfdScalar + Copy + SafeFromUsize + Sum> ResistanceAnalysis<T> {
    /// Create a new resistance analysis
    #[must_use]
    pub fn new() -> Self {
        Self {
            resistances: HashMap::new(),
            total_resistance: HydraulicResistance::from_base(T::ZERO),
            resistance_by_type: HashMap::new(),
            critical_paths: Vec::new(),
        }
    }

    /// Add resistance data for a component
    pub fn add_resistance(&mut self, id: String, resistance: HydraulicResistance<T>) {
        self.resistances.insert(id, resistance);
        self.total_resistance = HydraulicResistance::from_base(
            self.total_resistance.into_base() + resistance.into_base(),
        );
    }

    /// Add resistance by type
    pub fn add_resistance_by_type(
        &mut self,
        component_type: String,
        resistance: HydraulicResistance<T>,
    ) {
        let entry = self
            .resistance_by_type
            .entry(component_type)
            .or_insert(HydraulicResistance::from_base(T::ZERO));
        *entry = HydraulicResistance::from_base((*entry).into_base() + resistance.into_base());
    }

    /// Add a critical path
    pub fn add_critical_path(&mut self, path: Vec<String>) {
        self.critical_paths.push(path);
    }

    /// Get the average resistance
    pub fn average_resistance(&self) -> HydraulicResistance<T> {
        if self.resistances.is_empty() {
            HydraulicResistance::from_base(T::ZERO)
        } else {
            let sum: T = self.resistances.values().map(|r| r.into_base()).sum();
            HydraulicResistance::from_base(sum / T::from_usize_or_one(self.resistances.len()))
        }
    }

    /// Get the maximum resistance component
    pub fn max_resistance(&self) -> Option<(&String, HydraulicResistance<T>)> {
        self.resistances
            .iter()
            .max_by(|(_, a), (_, b)| {
                a.into_base()
                    .partial_cmp(&b.into_base())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(id, &r)| (id, r))
    }

    /// Get the minimum resistance component
    pub fn min_resistance(&self) -> Option<(&String, HydraulicResistance<T>)> {
        self.resistances
            .iter()
            .min_by(|(_, a), (_, b)| {
                a.into_base()
                    .partial_cmp(&b.into_base())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(id, &r)| (id, r))
    }

    /// Calculate parallel resistance
    pub fn parallel_resistance(&self) -> HydraulicResistance<T> {
        parallel_resistance(self.resistances.values().copied())
    }

    /// Calculate series resistance
    pub fn series_resistance(&self) -> HydraulicResistance<T> {
        series_resistance(self.resistances.values().copied())
    }
}

impl<T: CfdScalar + Copy + SafeFromUsize + Sum> Default for ResistanceAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}
