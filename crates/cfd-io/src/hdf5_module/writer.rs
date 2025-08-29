//! HDF5 writer for CFD datasets.

#[cfg(feature = "hdf5")]
use hdf5::{Dataset, File, Group, H5Type};

use crate::hdf5_module::metadata::DatasetMetadata;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use nalgebra::{DMatrix, DVector};
use std::path::Path;

/// HDF5 writer for large CFD datasets
pub struct Hdf5Writer {
    #[cfg(feature = "hdf5")]
    file: Option<File>,
    #[cfg(not(feature = "hdf5"))]
    _phantom: std::marker::PhantomData<()>,
}

impl Hdf5Writer {
    /// Create HDF5 writer instance
    pub fn create_instance() -> Self {
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
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
            ))
        }
    }

    /// Write vector dataset with chunking
    #[cfg(feature = "hdf5")]
    pub fn write_vector_chunked<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        data: &DVector<T>,
        chunk_size: usize,
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        let file = self
            .file
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

        let group = file
            .create_group(group_name)
            .or_else(|_| file.group(group_name))
            .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

        let dataset = group
            .current_dataset::<T>()
            .shape([data.len()])
            .chunk([chunk_size])
            .create(dataset_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to create dataset: {}", e)))?;

        // Write data in chunks
        data.as_slice()
            .chunks(chunk_size)
            .enumerate()
            .try_for_each(|(chunk_idx, chunk)| {
                let start_idx = chunk_idx * chunk_size;
                dataset
                    .write_slice(chunk, start_idx..)
                    .map_err(|e| Error::InvalidInput(format!("Failed to write chunk: {}", e)))
            })?;

        self.write_metadata_attributes(&dataset, metadata)?;
        Ok(())
    }

    #[cfg(not(feature = "hdf5"))]
    pub fn write_vector_chunked<T: RealField + Copy>(
        &self,
        _group_name: &str,
        _dataset_name: &str,
        _data: &DVector<T>,
        _chunk_size: usize,
        _metadata: &DatasetMetadata,
    ) -> Result<()> {
        Err(Error::InvalidConfiguration(
            "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
        ))
    }

    /// Write matrix dataset with 2D chunking
    #[cfg(feature = "hdf5")]
    pub fn write_matrix_chunked<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
        data: &DMatrix<T>,
        chunk_dims: (usize, usize),
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        let file = self
            .file
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

        let group = file
            .create_group(group_name)
            .or_else(|_| file.group(group_name))
            .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

        let dataset = group
            .current_dataset::<T>()
            .shape([data.nrows(), data.ncols()])
            .chunk([chunk_dims.0, chunk_dims.1])
            .create(dataset_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to create dataset: {}", e)))?;

        // Write matrix data
        let slice = data.as_slice();
        dataset
            .write_raw(slice)
            .map_err(|e| Error::InvalidInput(format!("Failed to write matrix data: {}", e)))?;

        self.write_metadata_attributes(&dataset, metadata)?;
        Ok(())
    }

    #[cfg(not(feature = "hdf5"))]
    pub fn write_matrix_chunked<T: RealField + Copy>(
        &self,
        _group_name: &str,
        _dataset_name: &str,
        _data: &DMatrix<T>,
        _chunk_dims: (usize, usize),
        _metadata: &DatasetMetadata,
    ) -> Result<()> {
        Err(Error::InvalidConfiguration(
            "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
        ))
    }

    /// Write metadata as HDF5 attributes
    #[cfg(feature = "hdf5")]
    fn write_metadata_attributes<T: H5Type>(
        &self,
        dataset: &Dataset,
        metadata: &DatasetMetadata,
    ) -> Result<()> {
        use hdf5::types::VarLenUnicode;

        dataset
            .current_attr::<VarLenUnicode>()
            .create("name")?
            .write_scalar(&metadata.name.clone().as_str().into())?;

        dataset
            .current_attr::<VarLenUnicode>()
            .create("data_type")?
            .write_scalar(&metadata.data_type.clone().as_str().into())?;

        if let Some(ref units) = metadata.units {
            dataset
                .current_attr::<VarLenUnicode>()
                .create("units")?
                .write_scalar(&units.clone().as_str().into())?;
        }

        if let Some(time_step) = metadata.time_step {
            dataset
                .current_attr::<f64>()
                .create("time_step")?
                .write_scalar(&time_step)?;
        }

        Ok(())
    }

    /// Close the HDF5 file
    pub fn close(&mut self) {
        #[cfg(feature = "hdf5")]
        {
            self.file = None;
        }
    }
}

impl Default for Hdf5Writer {
    fn default() -> Self {
        Self::create_instance()
    }
}
