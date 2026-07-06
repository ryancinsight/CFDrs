# cfd-io — Agent Reference

> **Role**: All file I/O for simulation data. No solver logic.  
> **Direct internal deps**: none. Solvers and validation crates depend on
> `cfd-io`; `cfd-io` does not depend back on solver, physics, mesh, or math
> crates.
> **Atlas providers**: Leto owns dense checkpoint/binary arrays, Eunomia owns
> scalar bounds, Consus owns HDF5, and RITK owns optional VTK.

---


<!-- AGENT-AUDIT-SNAPSHOT:START -->
## Verified Audit Snapshot (2026-07-04)

- Verified against `Cargo.toml`, `src/lib.rs`, and the top-level `src/` tree.
- Direct internal crate dependencies: none.
- Cargo features: `default`, `hdf5`, `vtk`.
- Provider dependencies: `leto`, `eunomia`, `consus-core`/`consus-hdf5`
  behind `hdf5`, and `ritk-vtk` behind `vtk`.
- `src/lib.rs` module surface: `binary`, `checkpoint`, `csv`, `error`,
  `hdf5` (feature), `vtk` (feature).
- Top-level `src/` entries: `binary.rs`, `checkpoint`, `csv`, `error.rs`,
  `formats.rs`, `hdf5/`, `json.rs`, `leto_arrays.rs`, `lib.rs`, `vtk.rs`.

<!-- AGENT-AUDIT-SNAPSHOT:END -->
## Purpose

`cfd-io` owns every read/write path for simulation data:

- **VTK**: structured and unstructured field output (`.vts`, `.vtu`, `.pvts`, `.pvtu`)  
- **CSV**: field data, convergence histories, validation tables  
- **Checkpoint**: binary snapshot with compression, metadata, integrity validation  
- **JSON**: metadata, configuration records  
- **Raw binary**: high-throughput field dumps  

No solver, mesh geometry, or physics logic lives in this crate.

---

## Module Structure

```
src/
  lib.rs                    pub mods; no prelude (hierarchical access preferred)
  error.rs                  Error, Result
  formats.rs                Format registry + auto-detection by extension
  json.rs                   Serde-based JSON read/write helpers
  leto_arrays.rs            row-major helpers for provider-owned Leto arrays
  binary.rs                 Raw binary Array1/Array2 payload dumps
  vtk.rs                    optional RITK VTK re-export

  csv/
    mod.rs
    types.rs                CsvRecord, ColumnSchema
    reader.rs               CsvReader — typed, header-aware
    writer.rs               CsvWriter — structured output with optional headers
    streaming.rs            StreamingReader/StreamingWriter — low-memory row-by-row export

  checkpoint/
    mod.rs
    data.rs                 Checkpoint — Leto-backed field snapshot container
    metadata.rs             CheckpointMetadata — timestamp, git hash, solver params
    manager.rs              CheckpointManager — save/load/rotate logic
    compression.rs          zstd / lz4 frame compression wrappers
    validator.rs            Checksum + schema validation on restore

  hdf5/
    mod.rs                  optional Consus-backed HDF5 writer
```

---

## Key APIs

### VTK Output

```rust
use cfd_io::vtk;

// `cfd_io::vtk::*` is the optional RITK VTK surface re-exported behind the
// `vtk` feature. Keep VTK model/writer behavior in RITK, not in cfd-io.
```

### Checkpoint

```rust
use cfd_io::checkpoint::{Checkpoint, CheckpointManager, CheckpointMetadata};
use leto::Array2;

let mut mgr = CheckpointManager::new("checkpoints/")?;
mgr.set_max_checkpoints(5);

let metadata = CheckpointMetadata::new(time, step, (ny, nx), (lx, ly));
let u = Array2::from_elem([ny, nx], 0.0);
let v = Array2::from_elem([ny, nx], 0.0);
let p = Array2::from_elem([ny, nx], 0.0);
let checkpoint = Checkpoint::new(metadata, u, v, p);

let path = mgr.save(&checkpoint)?;
let restored = mgr.load(path)?;
```

### CSV

```rust
use cfd_io::csv::{CsvWriter, StreamingWriter};

// Full table
CsvWriter::<f64>::new().write_time_series(
    "results/convergence.csv".as_ref(),
    &headers,
    rows,
)?;

// Streaming (append per time-step)
let mut stream =
    StreamingWriter::<f64>::create("results/residuals.csv".as_ref(), &["time", "residual"])?;
stream.write_row(&[time, residual])?;
```

---

## Design Rules

| Rule | Rationale |
|------|-----------|
| No `unwrap` on file paths | I/O paths are external inputs; always return `Result` |
| Streaming writes for large fields | Fields can exceed RAM when written all at once |
| Checkpoint validator runs on every load | Silent corruption is worse than a loud error |
| VTK XML preferred over legacy | Better tooling support in ParaView / VisIt |
| Parallel VTK uses reference files (`.pvts`, `.pvtu`) | Enables post-processing without re-running |
| Dense checkpoint and binary payloads use Leto arrays | Prevents nalgebra from re-entering the I/O boundary |
| Scalar bounds use Eunomia | Prevents direct `num-traits` drift in format code |

---

## Prohibited Patterns

- No physics computations — strip data from solver state before writing
- No in-memory caching of field data (fields are large; stream through)
- No direct dependency on `cfd-mesh`, `cfd-math`, `cfd-2d`, `cfd-3d` (one-way: solvers call cfd-io)
- No direct dependency on `cfd-core` only for error reporting; use the local
  file-format `Error`/`Result` type.
- No nalgebra or ndarray dense payloads; use Leto `Array1`/`Array2`.
- No direct `num-traits` scalar bounds; use Eunomia traits.



