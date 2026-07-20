# Atlas Crate Dependency Map

CFDrs consumes the Atlas stack through a layered dependency model:

- Numerics and traits: `eunomia`, `leto`, `leto-ops`
- SIMD and compute kernels: `hermes`, `hephaestus-wgpu`
- Runtime and orchestration: `moirai`
- Memory and allocation ecosystem: `mnemosyne`, `themis`, `melinoe`
- Spectral methods: `apollo-fft`

This appendix is the stable map for dependency-audit and migration work.

