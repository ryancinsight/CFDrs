# Summary

[Introduction](README.md)

---

# Part I — Foundations

- [CFDrs Architecture and Problem Setup](foundations.md)
  - [Example: cfd_demo](examples/cfd_demo.md)
  - [Example: enhanced_cfd_demo](examples/enhanced_cfd_demo.md)
  - [Example: adaptive_time_stepping_demo](examples/adaptive_time_stepping_demo.md)

# Part II — Core Flow Cases and Validation

- [Canonical Incompressible Benchmarks](core_flows.md)
  - [Example: cavity_validation](examples/cavity_validation.md)
  - [Example: pipe_flow_validation](examples/pipe_flow_validation.md)
  - [Example: turbulent_channel_flow](examples/turbulent_channel_flow.md)

# Part III — Turbulence and Multiphase

- [Turbulence Models and Cavitation](turbulence_multiphase.md)
  - [Example: turbulence_models_demo](examples/turbulence_models_demo.md)
  - [Example: turbulence_validation_demo](examples/turbulence_validation_demo.md)
  - [Example: venturi_cavitation](examples/venturi_cavitation.md)

# Part IV — Biomedical and Specialized Flows

- [Blood Flow and Rheology Workflows](biomedical_flows.md)
  - [Example: blood_flow_1d_validation](examples/blood_flow_1d_validation.md)
  - [Example: bifurcation_2d_blood_validation](examples/bifurcation_2d_blood_validation.md)
  - [Example: blood_rheology_models](examples/blood_rheology_models.md)

# Part V — Discretization and Solvers

- [Spectral, FEM, MUSCL, and Matrix-Free Methods](numerics_and_solvers.md)
  - [Example: spectral_3d_poisson](examples/spectral_3d_poisson.md)
  - [Example: fem_3d_stokes](examples/fem_3d_stokes.md)
  - [Example: matrix_free_demo](examples/matrix_free_demo.md)

# Part VI — Geometry, Meshing, and CSG

- [Geometry Construction and CFD Coupling](geometry_and_meshing.md)
  - [Example: csg_primitives_demo](examples/csg_primitives_demo.md)
  - [Example: csg_operations](examples/csg_operations.md)
  - [Example: csg_cfd_simulation](examples/csg_cfd_simulation.md)

# Part VII — Performance and Atlas Integration

- [SIMD, GPU, and Backend Migration](performance_and_atlas.md)
  - [Example: simd_performance_benchmark](examples/simd_performance_benchmark.md)
  - [Example: gpu_detection](examples/gpu_detection.md)
  - [Example: spectral_performance](examples/spectral_performance.md)

---

# Appendix

- [Atlas Crate Dependency Map](appendix_dependencies.md)
- [Migration Notes: ndarray/nalgebra/burn → leto/hephaestus/coeus](appendix_migration.md)
