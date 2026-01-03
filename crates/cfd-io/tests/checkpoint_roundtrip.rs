//! Bit-exact roundtrip tests for CheckpointManager

use cfd_io::checkpoint::{Checkpoint, CheckpointManager, CheckpointMetadata, CompressionStrategy};
use nalgebra::{dmatrix, DMatrix};
use proptest::prelude::*;
use tempfile::tempdir;

#[test]
fn checkpoint_roundtrip_no_compression() {
    let dir = tempdir().unwrap();
    let manager = CheckpointManager::new(dir.path()).unwrap();

    let ny = 10;
    let nx = 10;
    let metadata = CheckpointMetadata::new(1.234, 100, (ny, nx), (1.0, 1.0));

    let u = DMatrix::from_element(ny, nx, 1.0);
    let v = DMatrix::from_element(ny, nx, 2.0);
    let p = DMatrix::from_element(ny, nx, 3.0);

    let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

    let path = manager.save(&checkpoint).unwrap();
    let loaded = manager.load(&path).unwrap();

    assert_eq!(loaded.metadata.iteration, 100);
    assert_eq!(loaded.metadata.time.to_bits(), 1.234f64.to_bits());
    assert_eq!(loaded.u_velocity, u);
    assert_eq!(loaded.v_velocity, v);
    assert_eq!(loaded.pressure, p);
    assert_eq!(loaded.compute_checksum(), checkpoint.compute_checksum());
}

#[test]
fn checkpoint_roundtrip_zstd() {
    let dir = tempdir().unwrap();
    let mut manager = CheckpointManager::new(dir.path()).unwrap();
    manager.set_compression(CompressionStrategy::Zstd(3));

    let ny = 5;
    let nx = 5;
    let metadata = CheckpointMetadata::new(2.718, 200, (ny, nx), (2.0, 2.0));

    let u = dmatrix![1.0,2.0; 3.0,4.0; 5.0,6.0; 7.0,8.0; 9.0,10.0];
    let v = dmatrix![10.0,9.0; 8.0,7.0; 6.0,5.0; 4.0,3.0; 2.0,1.0];
    let p = DMatrix::from_element(5, 5, 42.0);

    let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

    let path = manager.save(&checkpoint).unwrap();
    let loaded = manager.load(&path).unwrap();

    assert_eq!(loaded.metadata.iteration, 200);
    assert_eq!(loaded.metadata.time.to_bits(), 2.718f64.to_bits());
    assert_eq!(loaded.u_velocity, u);
    assert_eq!(loaded.v_velocity, v);
    assert_eq!(loaded.pressure, p);
    assert_eq!(loaded.compute_checksum(), checkpoint.compute_checksum());
}

proptest::proptest! {
    #[test]
    fn proptest_roundtrip_small_grid(
        time in proptest::prelude::any::<f64>(),
        iteration in 0..1000usize,
        u_data in proptest::collection::vec(proptest::prelude::any::<f64>(), 1..25),
        v_data in proptest::collection::vec(proptest::prelude::any::<f64>(), 1..25),
        p_data in proptest::collection::vec(proptest::prelude::any::<f64>(), 1..25),
    ) {
        let n = (u_data.len() as f64).sqrt() as usize;
        if n * n != u_data.len() || n * n != v_data.len() || n * n != p_data.len() {
            return Ok(());
        }

        let ny = n;
        let nx = n;
        let metadata = CheckpointMetadata::new(
            time,
            iteration,
            (ny, nx),
            (1.0, 1.0),
        );

        let u = DMatrix::from_row_slice(ny, nx, &u_data);
        let v = DMatrix::from_row_slice(ny, nx, &v_data);
        let p = DMatrix::from_row_slice(ny, nx, &p_data);

        let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

        let dir = tempdir().unwrap();
        let manager = CheckpointManager::new(dir.path()).unwrap();

        let path = manager.save(&checkpoint).unwrap();
        let loaded = manager.load(&path).unwrap();

        prop_assert_eq!(loaded.metadata.iteration, iteration);
        prop_assert_eq!(loaded.metadata.time.to_bits(), time.to_bits());
        prop_assert_eq!(&loaded.u_velocity, &u);
        prop_assert_eq!(&loaded.v_velocity, &v);
        prop_assert_eq!(&loaded.pressure, &p);
        prop_assert_eq!(loaded.compute_checksum(), checkpoint.compute_checksum());
    }
}
