use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Mixture composition keyed by `fluid_id` with mass/volume fractions.
#[derive(Debug, Clone)]
pub struct MixtureComposition<T: RealField + Copy> {
    /// Fraction per fluid id.
    pub fractions: HashMap<i32, T>,
}

impl<T: RealField + Copy + FromPrimitive> MixtureComposition<T> {
    /// Create a new mixture and normalize to unit sum (when non-empty).
    #[must_use]
    pub fn new(mut fractions: HashMap<i32, T>) -> Self {
        let sum = fractions.values().copied().fold(T::zero(), |acc, v| acc + v);
        if sum > T::zero() {
            for value in fractions.values_mut() {
                *value /= sum;
            }
        }
        Self { fractions }
    }

    /// Empty mixture.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            fractions: HashMap::new(),
        }
    }

    /// Weighted blend of incoming mixtures.
    #[must_use]
    pub fn blend_weighted(inputs: &[(Self, T)]) -> Self {
        if inputs.is_empty() {
            return Self::empty();
        }

        let total_weight = inputs
            .iter()
            .map(|(_, w)| *w)
            .fold(T::zero(), |acc, v| acc + v);

        if total_weight <= T::zero() {
            return Self::empty();
        }

        let mut blended: HashMap<i32, T> = HashMap::new();
        for (mixture, weight) in inputs {
            for (fluid_id, frac) in &mixture.fractions {
                let contribution = (*frac) * (*weight) / total_weight;
                let entry = blended.entry(*fluid_id).or_insert(T::zero());
                *entry += contribution;
            }
        }

        Self::new(blended)
    }

    /// Compare with tolerance.
    #[must_use]
    pub fn approximately_equals(&self, other: &Self, tolerance: T) -> bool {
        let mut keys: Vec<i32> = self.fractions.keys().copied().collect();
        for key in other.fractions.keys() {
            if !keys.contains(key) {
                keys.push(*key);
            }
        }

        keys.into_iter().all(|k| {
            let a = *self.fractions.get(&k).unwrap_or(&T::zero());
            let b = *other.fractions.get(&k).unwrap_or(&T::zero());
            (a - b).abs() <= tolerance
        })
    }
}

/// Composition state at one simulation timepoint.
#[derive(Debug, Clone)]
pub struct CompositionState<T: RealField + Copy> {
    /// Simulation time.
    pub time: T,
    /// Node mixture compositions keyed by node index.
    pub node_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge mixture compositions keyed by edge index.
    pub edge_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge flow rates keyed by edge index for this snapshot.
    pub edge_flow_rates: HashMap<usize, T>,
}

impl<T: RealField + Copy> CompositionState<T> {
    /// Return average fluid concentrations in an edge at this state.
    ///
    /// In the current architecture, each edge stores a single mixed composition
    /// per snapshot, so this returns that composition as fluid concentrations.
    #[must_use]
    pub fn average_fluid_concentrations_in_edge(&self, edge_index: usize) -> Option<HashMap<i32, T>> {
        self.edge_mixtures
            .get(&edge_index)
            .map(|mixture| mixture.fractions.clone())
    }
}
