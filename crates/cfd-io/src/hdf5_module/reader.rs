//! HDF5 reader for CFD datasets.

#[cfg(feature = "hdf5")]
use hdf5::{File, H5Type};

use crate::hdf5_module::metadata::DatasetMetadata;
use cfd_core::{Error, Result};
use nalgebra::RealField;
use nalgebra::{DMatrix, DVector};
use std::path::Path;

/// HDF5 reader for large CFD datasets
pub struct Hdf5Reader {
    #[cfg(feature = "hdf5")]
    file: Option<File>,
    #[cfg(not(feature = "hdf5"))]
    _phantom: std::marker::PhantomData<()>,
}

impl Hdf5Reader {
    /// Create HDF5 reader instance
    pub fn create_instance() -> Self {
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
                "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
            ))
        }
    }

    /// Read vector dataset
    #[cfg(feature = "hdf5")]
    pub fn read_vector<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
    ) -> Result<DVector<T>> {
        let file = self
            .file
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

        let group = file
            .group(group_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

        let dataset = group
            .dataset(dataset_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access dataset: {}", e)))?;

        let shape = dataset.shape();
        if shape.len() != 1 {
            return Err(Error::InvalidInput(format!(
                "Expected 1D dataset, got {}D",
                shape.len()
            )));
        }

        let data = dataset
            .read_1d::<T>()
            .map_err(|e| Error::InvalidInput(format!("Failed to read data: {}", e)))?;

        Ok(DVector::from_vec(data.to_vec()))
    }

    #[cfg(not(feature = "hdf5"))]
    pub fn read_vector<T: RealField + Copy>(
        &self,
        _group_name: &str,
        _dataset_name: &str,
    ) -> Result<DVector<T>> {
        Err(Error::InvalidConfiguration(
            "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
        ))
    }

    /// Read matrix dataset
    #[cfg(feature = "hdf5")]
    pub fn read_matrix<T: RealField + Copy + H5Type>(
        &self,
        group_name: &str,
        dataset_name: &str,
    ) -> Result<DMatrix<T>> {
        let file = self
            .file
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

        let group = file
            .group(group_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

        let dataset = group
            .dataset(dataset_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access dataset: {}", e)))?;

        let shape = dataset.shape();
        if shape.len() != 2 {
            return Err(Error::InvalidInput(format!(
                "Expected 2D dataset, got {}D",
                shape.len()
            )));
        }

        let data = dataset
            .read_2d::<T>()
            .map_err(|e| Error::InvalidInput(format!("Failed to read data: {}", e)))?;

        let nrows = shape[0];
        let ncols = shape[1];
        let flat_data: Vec<T> = data.into_iter().flatten().collect();

        Ok(DMatrix::from_vec(nrows, ncols, flat_data))
    }

    #[cfg(not(feature = "hdf5"))]
    pub fn read_matrix<T: RealField + Copy>(
        &self,
        _group_name: &str,
        _dataset_name: &str,
    ) -> Result<DMatrix<T>> {
        Err(Error::InvalidConfiguration(
            "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
        ))
    }

    /// Read metadata from dataset
    #[cfg(feature = "hdf5")]
    pub fn read_metadata(&self, group_name: &str, dataset_name: &str) -> Result<DatasetMetadata> {
        use hdf5::types::VarLenUnicode;
        use std::collections::HashMap;

        let file = self
            .file
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("HDF5 file not opened".to_string()))?;

        let group = file
            .group(group_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access group: {}", e)))?;

        let dataset = group
            .dataset::<f64>(dataset_name)
            .map_err(|e| Error::InvalidInput(format!("Failed to access dataset: {}", e)))?;

        let name = dataset
            .attr("name")
            .ok()
            .and_then(|a| a.read_scalar::<VarLenUnicode>().ok())
            .map(|s| s.to_string())
            .unwrap_or_else(|| dataset_name.to_string());

        let data_type = dataset
            .attr("data_type")
            .ok()
            .and_then(|a| a.read_scalar::<VarLenUnicode>().ok())
            .map(|s| s.to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let units = dataset
            .attr("units")
            .ok()
            .and_then(|a| a.read_scalar::<VarLenUnicode>().ok())
            .map(|s| s.to_string());

        let time_step = dataset
            .attr("time_step")
            .ok()
            .and_then(|a| a.read_scalar::<f64>().ok());

        let dimensions = dataset.shape();

        Ok(DatasetMetadata {
            name,
            data_type,
            dimensions,
            units,
            time_step,
            attributes: HashMap::new(),
        })
    }

    #[cfg(not(feature = "hdf5"))]
    pub fn read_metadata(&self, _group_name: &str, _dataset_name: &str) -> Result<DatasetMetadata> {
        Err(Error::InvalidConfiguration(
            "HDF5 support not enabled. Enable 'hdf5' feature".to_string(),
        ))
    }

    /// Close the HDF5 file
    pub fn close(&mut self) {
        #[cfg(feature = "hdf5")]
        {
            self.file = None;
        }
    }
}

impl Default for Hdf5Reader {
    fn default() -> Self {
        Self::create_instance()
    }
}
