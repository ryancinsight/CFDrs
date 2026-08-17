use aequitas::systems::si::quantities::{Dimensionless, Time, VolumetricFlowRate};
use cfd_core::CfdScalar;
use std::collections::HashMap;

/// Canonical mixture key for transported RBC volume fraction.
pub const BLOOD_RBC_FLUID_ID: i32 = -10_001;
/// Canonical mixture key for transported plasma volume fraction.
pub const BLOOD_PLASMA_FLUID_ID: i32 = -10_002;

/// Mixture composition keyed by `fluid_id` with mass/volume fractions.
#[derive(Debug, Clone)]
pub struct MixtureComposition<T: CfdScalar + Copy> {
    /// Fraction per fluid id.
    pub fractions: HashMap<i32, Dimensionless<T>>,
}

impl<T: CfdScalar + Copy> MixtureComposition<T> {
    /// Create a new mixture and normalize to unit sum (when non-empty).
    #[must_use]
    pub fn new(mut fractions: HashMap<i32, Dimensionless<T>>) -> Self {
        let sum = fractions
            .values()
            .copied()
            .map(Dimensionless::into_base)
            .fold(T::ZERO, |acc, v| acc + v);
        if sum > T::ZERO {
            for value in fractions.values_mut() {
                *value = Dimensionless::from_base(value.into_base() / sum);
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

    /// Construct a two-species blood mixture from hematocrit.
    pub fn from_blood_hematocrit(hematocrit: Dimensionless<T>) -> Self {
        let hct = hematocrit.into_base().max(T::ZERO).min(T::ONE);
        let mut fractions = HashMap::new();
        fractions.insert(BLOOD_RBC_FLUID_ID, Dimensionless::from_base(hct));
        fractions.insert(
            BLOOD_PLASMA_FLUID_ID,
            Dimensionless::from_base(T::ONE - hct),
        );
        Self::new(fractions)
    }

    /// Return the transported blood hematocrit if the canonical RBC key exists.
    #[must_use]
    pub fn hematocrit(&self) -> Option<Dimensionless<T>> {
        self.fractions.get(&BLOOD_RBC_FLUID_ID).copied()
    }

    /// Weighted blend of incoming mixtures.
    #[must_use]
    pub fn blend_weighted(inputs: &[(&Self, Dimensionless<T>)]) -> Self {
        if inputs.is_empty() {
            return Self::empty();
        }

        let total_weight = inputs
            .iter()
            .map(|(_, w)| w.into_base())
            .fold(T::ZERO, |acc, v| acc + v);

        if total_weight <= T::ZERO {
            return Self::empty();
        }

        let mut blended: HashMap<i32, Dimensionless<T>> = HashMap::new();
        for (mixture, weight) in inputs {
            for (fluid_id, frac) in &mixture.fractions {
                let contribution = frac.into_base() * weight.into_base() / total_weight;
                let entry = blended
                    .entry(*fluid_id)
                    .or_insert_with(|| Dimensionless::from_base(T::ZERO));
                *entry = Dimensionless::from_base(entry.into_base() + contribution);
            }
        }

        Self::new(blended)
    }

    /// Weighted blend of owned incoming mixtures.
    pub(crate) fn blend_weighted_owned(inputs: &[(Self, T)]) -> Self {
        if inputs.is_empty() {
            return Self::empty();
        }

        let total_weight = inputs
            .iter()
            .map(|(_, w)| *w)
            .fold(T::ZERO, |acc, v| acc + v);

        if total_weight <= T::ZERO {
            return Self::empty();
        }

        let mut blended: HashMap<i32, Dimensionless<T>> = HashMap::new();
        for (mixture, weight) in inputs {
            for (fluid_id, frac) in &mixture.fractions {
                let contribution = frac.into_base() * *weight / total_weight;
                let entry = blended
                    .entry(*fluid_id)
                    .or_insert_with(|| Dimensionless::from_base(T::ZERO));
                *entry = Dimensionless::from_base(entry.into_base() + contribution);
            }
        }

        Self::new(blended)
    }

    /// Compare with tolerance.
    #[must_use]
    pub fn approximately_equals(&self, other: &Self, tolerance: Dimensionless<T>) -> bool {
        let mut keys: Vec<i32> = self.fractions.keys().copied().collect();
        for key in other.fractions.keys() {
            if !keys.contains(key) {
                keys.push(*key);
            }
        }

        keys.into_iter().all(|k| {
            let a = self
                .fractions
                .get(&k)
                .copied()
                .unwrap_or_else(|| Dimensionless::from_base(T::ZERO))
                .into_base();
            let b = other
                .fractions
                .get(&k)
                .copied()
                .unwrap_or_else(|| Dimensionless::from_base(T::ZERO))
                .into_base();
            <T as eunomia::NumericElement>::abs(a - b) <= tolerance.into_base()
        })
    }
}

/// Composition state at one simulation timepoint.
#[derive(Debug, Clone)]
pub struct CompositionState<T: CfdScalar + Copy> {
    /// Simulation time.
    pub time: Time<T>,
    /// Node mixture compositions keyed by node index.
    pub node_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge mixture compositions keyed by edge index.
    pub edge_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge flow rates keyed by edge index for this snapshot.
    pub edge_flow_rates: HashMap<usize, VolumetricFlowRate<T>>,
}

impl<T: CfdScalar + Copy> CompositionState<T> {
    /// Return average fluid concentrations in an edge at this state.
    ///
    /// In the current architecture, each edge stores a single mixed composition
    /// per snapshot, so this returns that composition as fluid concentrations.
    #[must_use]
    pub fn average_fluid_concentrations_in_edge(
        &self,
        edge_index: usize,
    ) -> Option<HashMap<i32, Dimensionless<T>>> {
        self.edge_mixtures
            .get(&edge_index)
            .map(|mixture| mixture.fractions.clone())
    }

    /// Return the transported hematocrit in a node mixture, if present.
    #[must_use]
    pub fn node_hematocrit(&self, node_index: usize) -> Option<Dimensionless<T>> {
        self.node_mixtures
            .get(&node_index)
            .and_then(MixtureComposition::hematocrit)
    }

    /// Return the transported hematocrit in an edge mixture, if present.
    #[must_use]
    pub fn edge_hematocrit(&self, edge_index: usize) -> Option<Dimensionless<T>> {
        self.edge_mixtures
            .get(&edge_index)
            .and_then(MixtureComposition::hematocrit)
    }
}
