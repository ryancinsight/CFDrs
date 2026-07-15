//! Bit-exact roundtrip tests for CheckpointManager

use cfd_io::checkpoint::{Checkpoint, CheckpointManager, CheckpointMetadata, CompressionStrategy};
use leto::Array2;
use proptest::prelude::*;
use std::f64::consts::E;
use tempfile::tempdir;

fn matrix_from_rows<const R: usize, const C: usize>(rows: [[f64; C]; R]) -> Array2<f64> {
    Array2::from_shape_vec([R, C], rows.into_iter().flatten().collect())
        .expect("test matrix shape and values must match")
}

fn assert_array_eq(left: &Array2<f64>, right: &Array2<f64>) {
    assert_eq!(left.shape(), right.shape());
    assert_eq!(left.clone().into_vec(), right.clone().into_vec());
}

#[test]
fn checkpoint_roundtrip_no_compression() {
    let dir = tempdir().unwrap();
    let manager = CheckpointManager::new(dir.path()).unwrap();

    let ny = 10;
    let nx = 10;
    let metadata = CheckpointMetadata::new(1.234, 100, (ny, nx), (1.0, 1.0));

    let u = Array2::from_elem([ny, nx], 1.0);
    let v = Array2::from_elem([ny, nx], 2.0);
    let p = Array2::from_elem([ny, nx], 3.0);

    let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

    let path = manager.save(&checkpoint).unwrap();
    let loaded = manager.load(&path).unwrap();

    assert_eq!(loaded.metadata.iteration, 100);
    assert_eq!(loaded.metadata.time.to_bits(), 1.234f64.to_bits());
    assert_array_eq(&loaded.u_velocity, &u);
    assert_array_eq(&loaded.v_velocity, &v);
    assert_array_eq(&loaded.pressure, &p);
    assert_eq!(loaded.compute_checksum(), checkpoint.compute_checksum());
}

#[test]
fn checkpoint_roundtrip_zstd() {
    let dir = tempdir().unwrap();
    let mut manager = CheckpointManager::new(dir.path()).unwrap();
    manager.set_compression(CompressionStrategy::Zstd(3));

    let ny = 5;
    let nx = 5;
    let metadata = CheckpointMetadata::new(E, 200, (ny, nx), (2.0, 2.0));

    let u = matrix_from_rows([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 10.0]]);
    let v = matrix_from_rows([[10.0, 9.0], [8.0, 7.0], [6.0, 5.0], [4.0, 3.0], [2.0, 1.0]]);
    let p = Array2::from_elem([5, 5], 42.0);

    let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

    let path = manager.save(&checkpoint).unwrap();
    let loaded = manager.load(&path).unwrap();

    assert_eq!(loaded.metadata.iteration, 200);
    assert_eq!(loaded.metadata.time.to_bits(), E.to_bits());
    assert_array_eq(&loaded.u_velocity, &u);
    assert_array_eq(&loaded.v_velocity, &v);
    assert_array_eq(&loaded.pressure, &p);
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

        let u = Array2::from_shape_vec([ny, nx], u_data.clone()).unwrap();
        let v = Array2::from_shape_vec([ny, nx], v_data.clone()).unwrap();
        let p = Array2::from_shape_vec([ny, nx], p_data.clone()).unwrap();

        let checkpoint = Checkpoint::new(metadata, u.clone(), v.clone(), p.clone());

        let dir = tempdir().unwrap();
        let manager = CheckpointManager::new(dir.path()).unwrap();

        let path = manager.save(&checkpoint).unwrap();
        let loaded = manager.load(&path).unwrap();

        prop_assert_eq!(loaded.metadata.iteration, iteration);
        prop_assert_eq!(loaded.metadata.time.to_bits(), time.to_bits());
        prop_assert_eq!(loaded.u_velocity.shape(), u.shape());
        prop_assert_eq!(loaded.v_velocity.shape(), v.shape());
        prop_assert_eq!(loaded.pressure.shape(), p.shape());
        prop_assert_eq!(loaded.u_velocity.clone().into_vec(), u.clone().into_vec());
        prop_assert_eq!(loaded.v_velocity.clone().into_vec(), v.clone().into_vec());
        prop_assert_eq!(loaded.pressure.clone().into_vec(), p.clone().into_vec());
        prop_assert_eq!(loaded.compute_checksum(), checkpoint.compute_checksum());
    }
}
