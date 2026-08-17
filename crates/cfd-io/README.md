# cfd-io

File I/O for the [CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

## What it provides

- **`csv`** — tabular field and time-series import/export.
- **`binary`** — compact binary field serialization.
- **`checkpoint`** — simulation state save and restore.
- **`hdf5`** (feature) — HDF5 output through the pure-Rust `consus-hdf5` backend.
- **`vtk`** (feature) — VTK read/write, re-exported from `ritk-vtk`.

Both optional backends are pure Rust. There is no MPI or collective-I/O layer;
output is single-process.

## Features

| Feature | Default | Effect |
|---|---|---|
| `hdf5` | no | HDF5 containers via `consus-core` + `consus-hdf5` |
| `vtk` | no | VTK creation and I/O via `ritk-vtk`, re-exported as `cfd_io::vtk` |

```bash
cargo build -p cfd-io --features hdf5,vtk
```

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
