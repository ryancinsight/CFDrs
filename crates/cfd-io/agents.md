# cfd-io — Agent Reference

> **Role**: All file I/O for simulation data. No solver logic.  
> **Depends on**: `cfd-core` only.

---

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
  error.rs                  IoError, IoResult
  formats.rs                Format registry + auto-detection by extension
  json.rs                   Serde-based JSON read/write helpers
  binary.rs                 Raw binary field dumps (f32/f64 array blobs)

  vtk/
    mod.rs
    types.rs                VtkDataSet, VtkArray, VtkCellType
    builder.rs              VtkBuilder — fluent API for structured/unstructured datasets
    writer.rs               VtkWriter — ASCII + binary VTK legacy + XML VTK
    reader.rs               VtkReader — legacy ASCII reader
    parallel.rs             Parallel VTK (PVD collection, PVTS/PVTU writers)

  csv/
    mod.rs
    types.rs                CsvRecord, ColumnSchema
    reader.rs               CsvReader — typed, header-aware
    writer.rs               CsvWriter — structured output with optional headers
    streaming.rs            StreamingCsvWriter — low-memory row-by-row export

  checkpoint/
    mod.rs
    data.rs                 CheckpointData — field snapshot container
    metadata.rs             CheckpointMetadata — timestamp, git hash, solver params
    manager.rs              CheckpointManager — save/load/rotate logic
    compression.rs          zstd / lz4 frame compression wrappers
    validator.rs            Checksum + schema validation on restore
```

---

## Key APIs

### VTK Output

```rust
use cfd_io::vtk::{VtkBuilder, VtkWriter};

// Structured grid (FDM/FVM output)
let dataset = VtkBuilder::structured_grid(nx, ny, nz)
    .with_points(coords)
    .add_cell_scalar("pressure", &p_field)
    .add_cell_vector("velocity", &u_field, &v_field, &w_field)
    .build()?;

VtkWriter::write_xml(&dataset, "output/step_0001.vts")?;

// Parallel output (MPI ranks write sub-domains)
use cfd_io::vtk::ParallelVtkWriter;
ParallelVtkWriter::write_piece(&dataset, rank, n_ranks, "output/step_0001")?;
```

### Checkpoint

```rust
use cfd_io::checkpoint::{CheckpointManager, CheckpointData};

let mut mgr = CheckpointManager::new("checkpoints/", max_keep: 5);

// Save
let data = CheckpointData::from_fields(&solver_state);
mgr.save(step, &data)?;

// Restore (last valid checkpoint)  
let restored = mgr.load_latest()?;
```

### CSV

```rust
use cfd_io::csv::{CsvWriter, StreamingCsvWriter};

// Full table
CsvWriter::write("results/convergence.csv", &headers, &rows)?;

// Streaming (append per time-step)
let mut stream = StreamingCsvWriter::open("results/residuals.csv", &["step","residual"])?;
stream.write_row(&[step.to_string(), residual.to_string()])?;
```

---

## Design Rules

| Rule | Rationale |
|------|-----------|
| No `unwrap` on file paths | I/O paths are external inputs; always return `IoResult` |
| Streaming writes for large fields | Fields can exceed RAM when written all at once |
| Checkpoint validator runs on every load | Silent corruption is worse than a loud error |
| VTK XML preferred over legacy | Better tooling support in ParaView / VisIt |
| Parallel VTK uses reference files (`.pvts`, `.pvtu`) | Enables post-processing without re-running |

---

## Prohibited Patterns

- No physics computations — strip data from solver state before writing
- No in-memory caching of field data (fields are large; stream through)
- No direct dependency on `cfd-mesh`, `cfd-math`, `cfd-2d`, `cfd-3d` (one-way: solvers call cfd-io)
