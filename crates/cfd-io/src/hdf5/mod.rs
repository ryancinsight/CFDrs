//! Real HDF5 output backed by the pure-Rust [`consus_hdf5`] implementation.
//!
//! This module writes standard, specification-conformant HDF5 files (verified
//! against `h5py` by the consus test-suite), in contrast to ad-hoc binary
//! containers. It is independent of MPI and usable on any platform.
//!
//! # Model
//!
//! A file is a flat root group of named, contiguous datasets plus a set of
//! root-group string attributes used for metadata. Datasets are written in the
//! native precision of the scalar type (`f32` → 32-bit, `f64` → 64-bit IEEE-754,
//! little-endian); no widening or narrowing is performed.
//!
//! # Example
//!
//! ```no_run
//! use cfd_io::hdf5::{write_hdf5, DatasetView};
//!
//! let pressure = [101_325.0_f64, 101_300.0, 101_280.0];
//! let datasets = [DatasetView { name: "pressure", shape: &[3], data: &pressure }];
//! let metadata = [("solver", "PISO"), ("units", "Pa")];
//! write_hdf5("field.h5", &datasets, &metadata).unwrap();
//! ```

use std::num::NonZeroUsize;
use std::path::Path;

use consus_core::{ByteOrder, Datatype, Shape, StringEncoding};
use consus_hdf5::file::writer::{DatasetCreationProps, FileCreationProps, Hdf5FileBuilder};

use crate::error::{Error, Result};

mod sealed {
    /// Prevents external implementations of [`super::Hdf5Float`], protecting the
    /// invariant that `BITS` matches the emitted little-endian byte width.
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

/// Scalar types writable as native-precision HDF5 floating-point datasets.
///
/// Sealed: only the IEEE-754 binary types this module can faithfully encode
/// implement it, so a dataset is always stored in the scalar's own precision.
pub trait Hdf5Float: Copy + sealed::Sealed {
    /// Bit width of the on-disk IEEE-754 representation (32 or 64).
    const BITS: usize;

    /// Append the little-endian byte encoding of `self` to `out`.
    fn push_le_bytes(self, out: &mut Vec<u8>);
}

impl Hdf5Float for f32 {
    const BITS: usize = 32;
    #[inline]
    fn push_le_bytes(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
}

impl Hdf5Float for f64 {
    const BITS: usize = 64;
    #[inline]
    fn push_le_bytes(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
}

/// A single named dataset to be written to the root group.
///
/// `shape` is the row-major extent of `data`; the product of its entries must
/// equal `data.len()`.
#[derive(Debug, Clone, Copy)]
pub struct DatasetView<'a, T: Hdf5Float> {
    /// Dataset name (root-group link name); must be non-empty.
    pub name: &'a str,
    /// Row-major dimension extents.
    pub shape: &'a [usize],
    /// Flat row-major values, length equal to the product of `shape`.
    pub data: &'a [T],
}

#[inline]
fn backend_err<E: core::fmt::Display>(e: E) -> Error {
    Error::Serialization(e.to_string())
}

/// Write `datasets` and `metadata` to a standard HDF5 file at `path`.
///
/// Each dataset is stored contiguously in the native precision of `T`. Each
/// metadata pair is stored as a UTF-8 fixed-length string attribute on the root
/// group.
///
/// # Errors
///
/// - [`Error::InvalidInput`] if a dataset name is empty, a dataset's declared
///   `shape` does not match its data length, or a metadata value is empty (a
///   zero-length HDF5 fixed string is not representable).
/// - [`Error::Serialization`] if the consus-hdf5 backend rejects an item.
/// - [`Error::Io`] if the file cannot be written.
pub fn write_hdf5<P, T>(
    path: P,
    datasets: &[DatasetView<'_, T>],
    metadata: &[(&str, &str)],
) -> Result<()>
where
    P: AsRef<Path>,
    T: Hdf5Float,
{
    // Native-precision IEEE-754 float datatype; `BITS` is a compile-time
    // constant of the scalar, so the cast to NonZeroUsize never fails.
    let bits = NonZeroUsize::new(T::BITS).expect("scalar bit width is non-zero");
    let datatype = Datatype::Float {
        bits,
        byte_order: ByteOrder::LittleEndian,
    };

    let mut builder = Hdf5FileBuilder::new(FileCreationProps::default());
    let dcpl = DatasetCreationProps::default();

    for ds in datasets {
        if ds.name.is_empty() {
            return Err(Error::InvalidInput("dataset name must be non-empty".into()));
        }
        let expected: usize = ds.shape.iter().product();
        if expected != ds.data.len() {
            return Err(Error::InvalidInput(format!(
                "dataset '{}': shape {:?} implies {} elements but data has {}",
                ds.name,
                ds.shape,
                expected,
                ds.data.len()
            )));
        }

        let mut raw = Vec::with_capacity(ds.data.len() * (T::BITS / 8));
        for &value in ds.data {
            value.push_le_bytes(&mut raw);
        }

        let shape = Shape::fixed(ds.shape);
        builder
            .add_dataset(ds.name, &datatype, &shape, &raw, &dcpl)
            .map_err(backend_err)?;
    }

    for (key, value) in metadata {
        let bytes = value.as_bytes();
        let length = NonZeroUsize::new(bytes.len()).ok_or_else(|| {
            Error::InvalidInput(format!("metadata '{key}': value must be non-empty"))
        })?;
        let attr_dt = Datatype::FixedString {
            length: length.get(),
            encoding: StringEncoding::Utf8,
        };
        builder
            .add_root_attribute(key, &attr_dt, &Shape::scalar(), bytes)
            .map_err(backend_err)?;
    }

    let image = builder.finish().map_err(backend_err)?;
    std::fs::write(path, image)?;
    Ok(())
}

#[cfg(all(test, feature = "hdf5"))]
mod tests {
    use super::*;
    use consus_hdf5::file::Hdf5File;
    use consus_io::MemCursor;

    /// Standard HDF5 superblock signature (\x89 H D F \r \n \x1a \n).
    const HDF5_SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n'];

    fn read_dataset_f64(file: &Hdf5File<MemCursor>, name: &str) -> Vec<f64> {
        let children = file.list_root_group().expect("list root group");
        let (_, addr, _) = children
            .iter()
            .find(|(n, _, _)| n == name)
            .unwrap_or_else(|| panic!("dataset '{name}' not found"));
        let dataset = file.dataset_at(*addr).expect("dataset header");
        let data_addr = dataset.data_address.expect("contiguous data address");
        let element_count: usize = dataset.shape.num_elements();
        let mut buf = vec![0u8; element_count * 8];
        file.read_contiguous_dataset_bytes(data_addr, 0, &mut buf)
            .expect("read dataset bytes");
        buf.chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }

    #[test]
    fn writes_a_genuine_hdf5_file_and_roundtrips_values_and_metadata() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("field.h5");

        // 2x3 row-major field plus a 1-D vector, distinct values per cell.
        let field: Vec<f64> = vec![1.5, -2.5, 3.0, 4.25, -5.0, 6.125];
        let velocity: Vec<f64> = vec![0.0, 10.0, -20.0, 30.5];
        let datasets = [
            DatasetView {
                name: "field",
                shape: &[2, 3],
                data: &field,
            },
            DatasetView {
                name: "velocity",
                shape: &[4],
                data: &velocity,
            },
        ];
        let metadata = [("solver", "PISO"), ("units", "SI")];

        write_hdf5(&path, &datasets, &metadata).expect("write hdf5");

        // The file must begin with the standard HDF5 signature — i.e. a real
        // HDF5 file, not an ad-hoc binary container.
        let bytes = std::fs::read(&path).expect("read back file");
        assert!(bytes.len() > 8, "file must be non-trivial");
        assert_eq!(&bytes[..8], &HDF5_SIGNATURE, "must be a real HDF5 file");

        let file = Hdf5File::open(MemCursor::from_bytes(bytes)).expect("open hdf5");

        // Exact value-semantic round-trip of both datasets.
        assert_eq!(read_dataset_f64(&file, "field"), field);
        assert_eq!(read_dataset_f64(&file, "velocity"), velocity);

        // Root-group metadata attributes round-trip exactly.
        let root = file.root_group();
        let attrs = file
            .attributes_at(root.object_header_address)
            .expect("read root attributes");
        for (key, value) in metadata {
            let attr = attrs
                .iter()
                .find(|a| a.name == key)
                .unwrap_or_else(|| panic!("metadata '{key}' missing"));
            assert_eq!(attr.raw_data, value.as_bytes(), "metadata '{key}' mismatch");
        }
    }

    #[test]
    fn native_precision_f32_datatype_is_32_bit() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("f32.h5");
        let data: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let datasets = [DatasetView {
            name: "x",
            shape: &[4],
            data: &data,
        }];
        write_hdf5(&path, &datasets, &[]).expect("write f32 hdf5");

        let bytes = std::fs::read(&path).expect("read back");
        let file = Hdf5File::open(MemCursor::from_bytes(bytes)).expect("open");
        let children = file.list_root_group().expect("list root");
        let (_, addr, _) = children
            .iter()
            .find(|(n, _, _)| n == "x")
            .expect("dataset x");
        let dataset = file.dataset_at(*addr).expect("dataset header");
        // Stored in native f32 precision (4 bytes/element), not widened to f64.
        match dataset.datatype {
            Datatype::Float { bits, .. } => assert_eq!(bits.get(), 32),
            other => panic!("expected Float datatype, got {other:?}"),
        }
    }

    #[test]
    fn shape_data_mismatch_is_rejected() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("bad.h5");
        let data: Vec<f64> = vec![1.0, 2.0, 3.0];
        let datasets = [DatasetView {
            name: "x",
            shape: &[2, 2],
            data: &data,
        }];
        let err = write_hdf5(&path, &datasets, &[]).unwrap_err();
        assert!(matches!(err, Error::InvalidInput(_)));
    }
}
