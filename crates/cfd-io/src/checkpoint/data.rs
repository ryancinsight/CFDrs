//! Checkpoint data structures

use crate::checkpoint::metadata::CheckpointMetadata;
use crate::leto_arrays::{row_major_values, try_for_each_row_major};
use eunomia::RealField;
use leto::Array2;
use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

/// Checkpoint data containing simulation state
#[derive(Debug, Clone)]
pub struct Checkpoint<T: RealField> {
    /// Metadata
    pub metadata: CheckpointMetadata,
    /// Velocity field (u component)
    pub u_velocity: Array2<T>,
    /// Velocity field (v component)
    pub v_velocity: Array2<T>,
    /// Pressure field
    pub pressure: Array2<T>,
    /// Temperature field (optional)
    pub temperature: Option<Array2<T>>,
    /// Turbulence kinetic energy (optional)
    pub turbulence_k: Option<Array2<T>>,
    /// Turbulence dissipation rate (optional)
    pub turbulence_epsilon: Option<Array2<T>>,
}

impl<T: RealField> Checkpoint<T> {
    /// Create a new checkpoint
    pub fn new(
        metadata: CheckpointMetadata,
        u_velocity: Array2<T>,
        v_velocity: Array2<T>,
        pressure: Array2<T>,
    ) -> Self {
        Self {
            metadata,
            u_velocity,
            v_velocity,
            pressure,
            temperature: None,
            turbulence_k: None,
            turbulence_epsilon: None,
        }
    }

    /// Validate checkpoint data consistency
    pub fn validate(&self) -> Result<(), String> {
        // Check dimensions consistency
        let (nx, ny) = self.metadata.dimensions;
        let [u_rows, u_cols] = self.u_velocity.shape();
        if u_cols != nx || u_rows != ny {
            return Err(format!(
                "U velocity dimensions mismatch: expected {ny}x{nx}, got {u_rows}x{u_cols}"
            ));
        }
        let [v_rows, v_cols] = self.v_velocity.shape();
        if v_cols != nx || v_rows != ny {
            return Err(format!(
                "V velocity dimensions mismatch: expected {ny}x{nx}, got {v_rows}x{v_cols}"
            ));
        }
        let [p_rows, p_cols] = self.pressure.shape();
        if p_cols != nx || p_rows != ny {
            return Err(format!(
                "Pressure dimensions mismatch: expected {ny}x{nx}, got {p_rows}x{p_cols}"
            ));
        }

        Ok(())
    }

    /// Get grid dimensions as (nx, ny)
    pub fn dimensions(&self) -> (usize, usize) {
        self.metadata.dimensions
    }
}

#[derive(Serialize, Deserialize)]
struct MatrixPayload<T> {
    shape: [usize; 2],
    values: Vec<T>,
}

impl<T: RealField> MatrixPayload<T> {
    fn from_array(array: &Array2<T>) -> Self {
        Self {
            shape: array.shape(),
            values: row_major_values(array),
        }
    }

    fn into_array(self) -> Result<Array2<T>, String> {
        Array2::from_shape_vec(self.shape, self.values)
            .map_err(|error| format!("invalid Leto checkpoint field payload: {error}"))
    }
}

#[derive(Serialize, Deserialize)]
struct CheckpointPayload<T> {
    metadata: CheckpointMetadata,
    u_velocity: MatrixPayload<T>,
    v_velocity: MatrixPayload<T>,
    pressure: MatrixPayload<T>,
    temperature: Option<MatrixPayload<T>>,
    turbulence_k: Option<MatrixPayload<T>>,
    turbulence_epsilon: Option<MatrixPayload<T>>,
}

impl<T> Serialize for Checkpoint<T>
where
    T: RealField + Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        CheckpointPayload {
            metadata: self.metadata.clone(),
            u_velocity: MatrixPayload::from_array(&self.u_velocity),
            v_velocity: MatrixPayload::from_array(&self.v_velocity),
            pressure: MatrixPayload::from_array(&self.pressure),
            temperature: self.temperature.as_ref().map(MatrixPayload::from_array),
            turbulence_k: self.turbulence_k.as_ref().map(MatrixPayload::from_array),
            turbulence_epsilon: self
                .turbulence_epsilon
                .as_ref()
                .map(MatrixPayload::from_array),
        }
        .serialize(serializer)
    }
}

impl<'de, T> Deserialize<'de> for Checkpoint<T>
where
    T: RealField + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let payload = CheckpointPayload::<T>::deserialize(deserializer)?;
        Ok(Self {
            metadata: payload.metadata,
            u_velocity: payload.u_velocity.into_array().map_err(D::Error::custom)?,
            v_velocity: payload.v_velocity.into_array().map_err(D::Error::custom)?,
            pressure: payload.pressure.into_array().map_err(D::Error::custom)?,
            temperature: payload
                .temperature
                .map(MatrixPayload::into_array)
                .transpose()
                .map_err(D::Error::custom)?,
            turbulence_k: payload
                .turbulence_k
                .map(MatrixPayload::into_array)
                .transpose()
                .map_err(D::Error::custom)?,
            turbulence_epsilon: payload
                .turbulence_epsilon
                .map(MatrixPayload::into_array)
                .transpose()
                .map_err(D::Error::custom)?,
        })
    }
}

impl Checkpoint<f64> {
    /// Computes deterministic u128 checksum for bit-exact verification.
    ///
    /// # Invariants
    /// * Roundtrip preservation: `checkpoint.compute_checksum() == loaded.compute_checksum()`
    /// * Parallel consistency: `allreduce_xor(local_checksums) == global_checksum`
    ///
    /// Uses dual DefaultHasher on metadata fields (skip checksum) + row-major matrix bits.
    pub fn compute_checksum(&self) -> u128 {
        let mut hasher1 = DefaultHasher::new();
        let mut hasher2 = DefaultHasher::new();

        let meta = &self.metadata;
        hasher1.write_u32(meta.version);
        hasher1.write_u64(meta.time.to_bits());
        hasher1.write_u64(meta.iteration as u64);
        hasher1.write_u64(meta.dimensions.0 as u64);
        hasher1.write_u64(meta.dimensions.1 as u64);
        hasher1.write_u64(meta.domain_size.0.to_bits());
        hasher1.write_u64(meta.domain_size.1.to_bits());
        hasher1.write_u64(meta.timestamp);

        hasher2.write_u32(meta.version);
        hasher2.write_u64(meta.time.to_bits());
        hasher2.write_u64(meta.iteration as u64);
        hasher2.write_u64(meta.dimensions.0 as u64);
        hasher2.write_u64(meta.dimensions.1 as u64);
        hasher2.write_u64(meta.domain_size.0.to_bits());
        hasher2.write_u64(meta.domain_size.1.to_bits());
        hasher2.write_u64(meta.timestamp);

        fn hash_matrix(
            matrix: &Array2<f64>,
            hasher1: &mut DefaultHasher,
            hasher2: &mut DefaultHasher,
        ) {
            try_for_each_row_major(matrix, |value| {
                let bits = value.to_bits();
                hasher1.write_u64(bits);
                hasher2.write_u64(bits);
                Ok::<(), core::convert::Infallible>(())
            })
            .expect("invariant: infallible checksum traversal cannot fail");
        }

        hash_matrix(&self.u_velocity, &mut hasher1, &mut hasher2);
        hash_matrix(&self.v_velocity, &mut hasher1, &mut hasher2);
        hash_matrix(&self.pressure, &mut hasher1, &mut hasher2);
        if let Some(ref mat) = self.temperature {
            hash_matrix(mat, &mut hasher1, &mut hasher2);
        }
        if let Some(ref mat) = self.turbulence_k {
            hash_matrix(mat, &mut hasher1, &mut hasher2);
        }
        if let Some(ref mat) = self.turbulence_epsilon {
            hash_matrix(mat, &mut hasher1, &mut hasher2);
        }

        let low = hasher1.finish();
        let high = hasher2.finish();
        u128::from(low) | (u128::from(high) << 64)
    }
}
