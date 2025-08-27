//! Data chunking for streaming I/O operations.

use crate::hdf5_module::metadata::DatasetMetadata;
use nalgebra::RealField;

/// Large dataset chunk for streaming I/O
#[derive(Debug, Clone)]
pub struct DataChunk<T: RealField + Copy> {
    /// Chunk data
    pub data: Vec<T>,
    /// Chunk offset in global dataset
    pub offset: Vec<usize>,
    /// Chunk dimensions
    pub dimensions: Vec<usize>,
    /// Metadata
    pub metadata: DatasetMetadata,
}

impl<T: RealField + Copy> DataChunk<T> {
    /// Create a data chunk
    pub fn create(
        data: Vec<T>,
        offset: Vec<usize>,
        dimensions: Vec<usize>,
        metadata: DatasetMetadata,
    ) -> Self {
        Self {
            data,
            offset,
            dimensions,
            metadata,
        }
    }

    /// Get the total size of the chunk
    pub fn size(&self) -> usize {
        self.dimensions.iter().product()
    }

    /// Validate that data size matches dimensions
    pub fn validate(&self) -> bool {
        self.data.len() == self.size()
    }

    /// Get a slice of the data
    pub fn data_slice(&self) -> &[T] {
        &self.data
    }

    /// Get mutable slice of the data
    pub fn data_slice_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}

/// Chunking strategy for large datasets
#[derive(Debug, Clone)]
pub struct ChunkingStrategy {
    /// Maximum chunk size in elements
    pub max_chunk_size: usize,
    /// Preferred chunk dimensions
    pub preferred_dimensions: Option<Vec<usize>>,
}

impl Default for ChunkingStrategy {
    fn default() -> Self {
        Self {
            max_chunk_size: 1_000_000, // 1M elements default
            preferred_dimensions: None,
        }
    }
}

impl ChunkingStrategy {
    /// Create a strategy with specified max chunk size
    pub fn with_max_size(max_chunk_size: usize) -> Self {
        Self {
            max_chunk_size,
            preferred_dimensions: None,
        }
    }

    /// Set preferred dimensions for chunks
    pub fn with_dimensions(mut self, dimensions: Vec<usize>) -> Self {
        self.preferred_dimensions = Some(dimensions);
        self
    }

    /// Calculate optimal chunk dimensions for a dataset
    pub fn calculate_chunk_dimensions(&self, dataset_dimensions: &[usize]) -> Vec<usize> {
        if let Some(ref preferred) = self.preferred_dimensions {
            return preferred.clone();
        }

        // Calculate balanced chunk dimensions
        let ndims = dataset_dimensions.len();
        let chunk_size_per_dim = (self.max_chunk_size as f64).powf(1.0 / ndims as f64) as usize;

        dataset_dimensions
            .iter()
            .map(|&dim| dim.min(chunk_size_per_dim))
            .collect()
    }
}
