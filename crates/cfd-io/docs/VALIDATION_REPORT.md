# cfd-io Validation Report

## Modules
- HDF5 reader/writer with chunking/compression; VTK I/O (reader/writer/types); CSV/JSON; checkpoint manager and validator.

## Tests
- Planned round-trip serialization tests (f32/f64) for field arrays and metadata; schema validation for boundary and mesh metadata; streaming correctness.

## Invariants
- Dimensional fields preserved across formats; metadata integrity; deterministic serialization (within format constraints).

## Units
- SI units for physical fields; documented conversions when applicable.

## Assumptions
- File format compatibility; compression settings and chunking aligned with performance and fidelity goals.

## References
- VTK/HDF5 format specifications; best practices for scientific data I/O.
