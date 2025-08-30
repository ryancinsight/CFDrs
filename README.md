# CFD Suite - Rust Implementation

**Version 0.1.0-ALPHA** - Not Production Ready

## ⚠️ HONEST STATUS WARNING ⚠️

This codebase is in **ALPHA** state. While the architecture is sound and core physics implementations are validated, significant work remains before production use.

### What Works ✅
- Core library compiles without errors
- Physics implementations validated against literature:
  - k-ε turbulence model (Launder & Spalding 1974)
  - k-ω SST model (Menter 1994)
  - Vorticity-stream formulation (Anderson 1995)
- Zero-copy patterns in critical paths
- Modular architecture (<500 lines per module)
- Basic SIMD infrastructure (runtime detection)

### What's Broken ❌
- **Many examples don't compile** - API has evolved, examples haven't
- **GPU support is feature-gated** - Not enabled by default
- **SIMD not integrated** - Infrastructure exists but unused in solvers
- **Limited test coverage** - Tests compile but coverage is minimal
- **Incomplete documentation** - 27+ missing documentation warnings

### What's Missing 🔴
- Performance benchmarks
- Comprehensive test suite
- User documentation
- API stability guarantees
- Production validation

## Architecture

8 specialized crates in a workspace:
- `cfd-core` - Core abstractions and traits
- `cfd-math` - Numerical methods and linear algebra
- `cfd-mesh` - Mesh generation and manipulation
- `cfd-io` - File I/O (VTK, HDF5 when enabled)
- `cfd-1d` - 1D solvers (pipe networks, channels)
- `cfd-2d` - 2D solvers (SIMPLE, PISO, vorticity-stream)
- `cfd-3d` - 3D solvers (FEM Stokes, spectral methods)
- `cfd-validation` - Validation and verification tools

## Building

```bash
# Basic build
cargo build --workspace

# With GPU support (experimental)
cargo build --workspace --features gpu

# Run tests (limited coverage)
cargo test --workspace

# Build documentation
cargo doc --workspace --open
```

## Examples

⚠️ **Warning**: Many examples are currently broken due to API changes.

Working examples:
- `cargo run --example cfd_demo`
- `cargo run --example 2d_heat_diffusion`

## Performance

### Claimed vs Reality
- **Claimed**: "Zero-copy throughout" → **Reality**: Mostly achieved in hot paths
- **Claimed**: "SIMD vectorization" → **Reality**: Infrastructure only, not integrated
- **Claimed**: "GPU acceleration" → **Reality**: Feature-gated, untested
- **Claimed**: "100k cells/sec" → **Reality**: No benchmarks to verify

## Contributing

This codebase needs:
1. **Fix broken examples** - Update to current API
2. **Complete documentation** - Document all public APIs
3. **Add benchmarks** - Validate performance claims
4. **Increase test coverage** - Currently minimal
5. **Integrate SIMD** - Use the existing infrastructure
6. **Enable GPU by default** - If stable enough

## Physics Validation

The following implementations have been validated against literature:

| Component | Reference | Status |
|-----------|-----------|--------|
| k-ε model | Launder & Spalding (1974) | ✅ Validated |
| k-ω SST | Menter (1994) | ✅ Validated |
| Vorticity-Stream | Anderson (1995) | ✅ Validated |
| SIMPLE algorithm | Patankar (1980) | ⚠️ Partial |
| VOF method | Hirt & Nichols (1981) | ⚠️ Untested |

## License

MIT OR Apache-2.0

## Acknowledgments

This is an honest attempt at a CFD suite in Rust. It's not production-ready, but the foundation is solid. The architecture is clean, the physics is mostly correct, and the performance patterns are sound. It just needs more work to be actually usable.

**Bottom Line**: Use this for learning and experimentation, not production CFD work. Yet.