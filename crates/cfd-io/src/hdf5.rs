//! HDF5 support for large dataset I/O operations.

#[cfg(feature = "hdf5")]
use hdf5::{File, Group, Dataset, H5Type};
use cfd_core::{Error, Result, RealField};
use nalgebra::{DVector, DMatrix};
use std::path::Path;
use serde::{Serialize, Deserialize};

/// HDF5 writer for large CFD datasets
pub struct Hdf5Writer {
    #[cfg(feature = "hdf5")]
    file: Option<File>,
    #[cfg(not(feature = "hdf5"))]
    _phantom: std::marker::PhantomData<()>,
}

/// HDF5 reader for large CFD datasets
pub struct Hdf5Reader {
    #[cfg(feature = "hdf5")]
    file: Option<File>,
    #[cfg(not(feature = "hdf5"))]
    _phantom: std::marker::PhantomData<()>,
}

/// Metadata for CFD datasets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetMetadata {
    /// Dataset name
    pub name: String,
    /// Data type description
    pub data_type: String,
    /// Dimensions
    pub dimensions: Vec<usize>,
    /// Physical units
    pub units: Option<String>,
    /// Time step (if time-dependent)
    pub time_step: Option<f64>,
    /// Additional attributes
    pub attributes: std::collections::HashMap<String, String>,
}

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

impl Hdf5Writer {
    /// Create new HDF5 writer
    pub fn new() -> Self {
        Self {
            #[cfg(feature = "hdf5")]
            file: None,
            #[cfg(not(feature = "hdf5"))]
            _phantom: std::marker::PhantomData,
        }
    }

    /// Open HDF5 file for writing
    pub fn open<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            let file = File::create(path.as_ref())
                .map_err(|e| Error::InvalidInput(format!("Failed to create HDF5 file: {}", e)))?;
            self.file = Some(file);
            Ok(())
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Write large vector dataset with chunking
    pub fn write_vector_chunked<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        data: &DVector<T>,
        chunk_size: usize,
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            let file = self.file.as_ref()
                .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

            // Create group if it doesn't exist
            let group = file.create_group(group_name)
                .or_else(|_| file.group(group_name))
                .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

            // Create dataset with chunking
            let dataset = group
                .current_dataset::<T>()
                .shape([data.len()])
                .chunk([chunk_size])
                .create(dataset_name)
                .map_err(|e| Error::InvalidInput(format!("Failed to create dataset: {}", e)))?;

            // Write data in chunks using iterator
            data.as_slice()
                .chunks(chunk_size)
                .enumerate()
                .try_for_each(|(chunk_idx, chunk)| {
                    let start_idx = chunk_idx * chunk_size;
                    dataset.write_slice(chunk, start_idx..)
                        .map_err(|e| Error::InvalidInput(format!("Failed to write chunk: {}", e)))
                })?;

            // Write metadata as attributes
            self.write_metadata_attributes(&dataset, metadata)?;

            Ok(())
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Write large matrix dataset with 2D chunking
    pub fn write_matrix_chunked<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        data: &DMatrix<T>,
        chunk_dims: (usize, usize),
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            let file = self.file.as_ref()
                .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

            let group = file.create_group(group_name)
                .or_else(|_| file.group(group_name))
                .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

            // Create 2D dataset with chunking
            let dataset = group
                .current_dataset::<T>()
                .shape([data.nrows(), data.ncols()])
                .chunk([chunk_dims.0, chunk_dims.1])
                .create(dataset_name)
                .map_err(|e| Error::InvalidInput(format!("Failed to create dataset: {}", e)))?;

            // Write matrix data
            let matrix_data: Vec<T> = data.iter().cloned().collect();
            dataset.write(&matrix_data)
                .map_err(|e| Error::InvalidInput(format!("Failed to write matrix: {}", e)))?;

            // Write metadata
            self.write_metadata_attributes(&dataset, metadata)?;

            Ok(())
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Write time series data with compression
    pub fn write_time_series<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        time_steps: &[f64],
        data: &[DVector<T>],
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            if time_steps.len() != data.len() {
                return Err(Error::InvalidConfiguration(
                    "Time steps and data length mismatch".to_string()
                ));
            }

            let file = self.file.as_ref()
                .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

            let group = file.create_group(group_name)
                .or_else(|_| file.group(group_name))
                .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

            // Write time steps
            let time_dataset = group
                .current_dataset::<f64>()
                .shape([time_steps.len()])
                .create("time_steps")
                .map_err(|e| Error::InvalidInput(format!("Failed to create time dataset: {}", e)))?;
            
            time_dataset.write(time_steps)
                .map_err(|e| Error::InvalidInput(format!("Failed to write time steps: {}", e)))?;

            // Write data for each time step
            let data_len = data.first().map(|v| v.len()).unwrap_or(0);
            let flattened_data: Vec<T> = data
                .iter()
                .flat_map(|vec| vec.iter().cloned())
                .collect();

            let data_dataset = group
                .current_dataset::<T>()
                .shape([data.len(), data_len])
                .chunk([1, data_len])
                .deflate(6) // Compression level
                .create(dataset_name)
                .map_err(|e| Error::InvalidInput(format!("Failed to create data dataset: {}", e)))?;

            data_dataset.write(&flattened_data)
                .map_err(|e| Error::InvalidInput(format!("Failed to write time series data: {}", e)))?;

            // Write metadata
            self.write_metadata_attributes(&data_dataset, metadata)?;

            Ok(())
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Write metadata as HDF5 attributes
    #[cfg(feature = "hdf5")]
    fn write_metadata_attributes<T: H5Type>(
        &self,
        dataset: &Dataset,
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        // Write basic metadata
        dataset.current_attr::<hdf5::types::VarLenUnicode>()
            .create("name")?
            .write_scalar(&metadata.name.clone().as_str().into())?;

        dataset.current_attr::<hdf5::types::VarLenUnicode>()
            .create("data_type")?
            .write_scalar(&metadata.data_type.as_str().into())?;

        // Write dimensions
        let dims: Vec<u64> = metadata.dimensions.iter().map(|&d| d as u64).collect();
        dataset.current_attr::<u64>()
            .shape(&[dims.len()])
            .create("dimensions")?
            .write(&dims)?;

        // Write optional fields
        if let Some(ref units) = metadata.units {
            dataset.current_attr::<hdf5::types::VarLenUnicode>()
                .create("units")?
                .write_scalar(&units.as_str().into())?;
        }

        if let Some(time_step) = metadata.time_step {
            dataset.current_attr::<f64>()
                .create("time_step")?
                .write_scalar(&time_step)?;
        }

        // Write additional attributes
        for (key, value) in &metadata.attributes {
            dataset.current_attr::<hdf5::types::VarLenUnicode>()
                .create(key)?
                .write_scalar(&value.as_str().into())?;
        }

        Ok(())
    }

    /// Close the HDF5 file
    pub fn close(&mut self) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            self.file = None;
        }
        Ok(())
    }
}

impl Hdf5Reader {
    /// Create new HDF5 reader
    pub fn new() -> Self {
        Self {
            #[cfg(feature = "hdf5")]
            file: None,
            #[cfg(not(feature = "hdf5"))]
            _phantom: std::marker::PhantomData,
        }
    }

    /// Open HDF5 file for reading
    pub fn open<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            let file = File::open(path.as_ref())
                .map_err(|e| Error::InvalidInput(format!("Failed to open HDF5 file: {}", e)))?;
            self.file = Some(file);
            Ok(())
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Read vector dataset in chunks
    pub fn read_vector_chunked<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        chunk_size: usize,
    ) -> Result<impl Iterator<Item = Result<DataChunk<T>>>> {
        #[cfg(feature = "hdf5")]
        {
            let file = self.file.as_ref()
                .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

            let group = file.group(group_name)
                .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

            let dataset = group.dataset(dataset_name)
                .map_err(|e| Error::InvalidInput(format!("Failed to access dataset: {}", e)))?;

            let shape = dataset.shape();
            let total_size = shape[0];
            
            // Read metadata
            let metadata = self.read_metadata_attributes(&dataset)?;

            // Create chunk iterator
            let chunk_iter = (0..total_size)
                .step_by(chunk_size)
                .map(move |start_idx| {
                    let end_idx = (start_idx + chunk_size).min(total_size);
                    let chunk_data: Vec<T> = dataset
                        .read_slice_1d(start_idx..end_idx)
                        .map_err(|e| Error::InvalidInput(format!("Failed to read chunk: {}", e)))?;

                    Ok(DataChunk {
                        data: chunk_data,
                        offset: vec![start_idx],
                        dimensions: vec![end_idx - start_idx],
                        metadata: metadata,
                    })
                });

            Ok(chunk_iter)
        }
        #[cfg(not(feature = "hdf5"))]
        {
            Err(Error::InvalidConfiguration(
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string()
            ))
        }
    }

    /// Read metadata from dataset attributes with proper error handling
    #[cfg(feature = "hdf5")]
    fn read_metadata_attributes(&self, dataset: &Dataset) -> Result<DatasetMetadata> {
        // Read required attributes with proper error handling
        let name = dataset.attr("name")
            .and_then(|attr| attr.read_scalar::<hdf5::types::VarLenUnicode>())
            .map(|s| s.to_string())
            .map_err(|e| Error::InvalidInput(format!("Failed to read required 'name' attribute: {}", e)))?;

        let data_type = dataset.attr("data_type")
            .and_then(|attr| attr.read_scalar::<hdf5::types::VarLenUnicode>())
            .map(|s| s.to_string())
            .map_err(|e| Error::InvalidInput(format!("Failed to read required 'data_type' attribute: {}", e)))?;

        let dimensions = dataset.attr("dimensions")
            .and_then(|attr| attr.read::<u64>())
            .map(|dims| dims.into_iter().map(|d| d as usize).collect())
            .map_err(|e| Error::InvalidInput(format!("Failed to read required 'dimensions' attribute: {}", e)))?;

        // Read optional attributes
        let units = dataset.attr("units")
            .and_then(|attr| attr.read_scalar::<hdf5::types::VarLenUnicode>())
            .map(|s| Some(s.to_string()))
            .unwrap_or(None);

        let time_step = dataset.attr("time_step")
            .and_then(|attr| attr.read_scalar::<f64>())
            .ok();

        // Read all additional attributes into the attributes map
        let mut attributes = std::collections::HashMap::new();

        // Get all attribute names from the dataset
        let attr_names = dataset.attr_names()
            .map_err(|e| Error::InvalidInput(format!("Failed to get attribute names: {}", e)))?;

        // Skip the standard attributes we've already processed
        let standard_attrs = ["name", "data_type", "dimensions", "units", "time_step"];

        for attr_name in attr_names {
            if !standard_attrs.contains(&attr_name.as_str()) {
                // Try to read as string attribute
                if let Ok(attr) = dataset.attr(&attr_name) {
                    if let Ok(value) = attr.read_scalar::<hdf5::types::VarLenUnicode>() {
                        attributes.insert(attr_name, value.to_string());
                    } else if let Ok(value) = attr.read_scalar::<f64>() {
                        // Convert numeric values to strings for consistency
                        attributes.insert(attr_name, value.to_string());
                    } else if let Ok(value) = attr.read_scalar::<i64>() {
                        attributes.insert(attr_name, value.to_string());
                    }
                    // Add more type conversions as needed
                }
            }
        }

        Ok(DatasetMetadata {
            name,
            data_type,
            dimensions,
            units,
            time_step,
            attributes,
        })
    }

    /// Close the HDF5 file
    pub fn close(&mut self) -> Result<()> {
        #[cfg(feature = "hdf5")]
        {
            self.file = None;
        }
        Ok(())
    }
}

impl Default for Hdf5Writer {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for Hdf5Reader {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;

    #[test]
    fn test_hdf5_writer_creation() {
        let writer = Hdf5Writer::new();
        assert!(writer.file.is_none());
    }

    #[test]
    fn test_hdf5_reader_creation() {
        let reader = Hdf5Reader::new();
        assert!(reader.file.is_none());
    }

    #[test]
    fn test_metadata_creation() {
        let metadata = DatasetMetadata {
            name: "test_dataset".to_string(),
            data_type: "f64".to_string(),
            dimensions: vec![100, 50],
            units: Some("m/s".to_string()),
            time_step: Some(0.01),
            attributes: std::collections::HashMap::new(),
        };

        assert_eq!(metadata.name.clone(), "test_dataset");
        assert_eq!(metadata.dimensions, vec![100, 50]);
    }
}
