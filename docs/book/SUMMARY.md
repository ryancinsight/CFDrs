# Summary

[Introduction](README.md)

# Part I — Governing Equations and Physical Models

- [1. CFDrs Architecture and Problem Setup](foundations.md)
  - [Example: CFD Demo](examples/cfd_demo.md)
  - [Example: Enhanced CFD Demo](examples/enhanced_cfd_demo.md)
- [2. Canonical Incompressible Benchmarks](core_flows.md)
- [3. Governing Equations and Discretization](governing_equations.md)
- [4. Pressure-Velocity Coupling and Time Integration](pressure_velocity.md)
  - [Example: SIMPLEC PIMPLE Demo](examples/simplec_pimple_demo.md)
  - [Example: Adaptive Time Stepping Demo](examples/adaptive_time_stepping_demo.md)

# Part II — Turbulence and Multiphase Physics

- [5. Turbulence Models: k-ε, k-ω SST, Smagorinsky LES](turbulence_multiphase.md)
  - [Example: Turbulence Models Demo](examples/turbulence_models_demo.md)
  - [Example: Turbulence Validation Demo](examples/turbulence_validation_demo.md)
  - [Example: Turbulence Momentum Integration Demo](examples/turbulence_momentum_integration_demo.md)
  - [Example: Turbulent Channel Flow](examples/turbulent_channel_flow.md)
- [6. Cavitation: Liquid-Vapour Phase Transition](cavitation.md)
  - [Example: Venturi Cavitation](examples/venturi_cavitation.md)
  - [Example: Simple Cavitation](examples/simple_cavitation.md)
  - [Example: Cavitation Damage Simulation](examples/cavitation_damage_simulation.md)
  - [Example: 1D Venturi Cavitation](examples/1d_venturi_cavitation.md)

# Part III — Biomedical and Rheological Flows

- [7. Non-Newtonian Fluids and Blood Rheology](biomedical_flows.md)
  - [Example: Blood Rheology Models](examples/blood_rheology_models.md)
  - [Example: Blood Flow 1D Validation](examples/blood_flow_1d_validation.md)
- [8. Vascular Bifurcations and Stenosis](vascular_bifurcations.md)
  - [Example: Bifurcation 2D Blood Validation](examples/bifurcation_2d_blood_validation.md)
  - [Example: Venturi Blood Flow Validation](examples/venturi_blood_flow_validation.md)
  - [Example: Cross Fidelity Branching](examples/cross_fidelity_branching.md)
- [9. Microfluidics and Millifluidic Networks](crate_1d_flows.md)
  - [Example: Blood Bifurcation](examples/blood_bifurcation.md)
  - [Example: Microfluidic Chip](examples/microfluidic_chip.md)
  - [Example: FDA Shear Limit Screening](examples/fda_shear_limit_screening.md)
  - [Example: TPMS Blood 1D](examples/tpms_blood_1d.md)
  - [Example: Medical Millifluidic Screening](examples/medical_millifluidic_screening.md)
  - [Example: Hemolysis Serpentine Analysis](examples/hemolysis_serpentine_analysis.md)
  - [Example: Serpentine Mixing Comprehensive](examples/serpentine_mixing_comprehensive.md)
  - [Example: Venturi Parallel Analysis](examples/venturi_parallel_analysis.md)

# Part IV — Numerical Methods and Linear Solvers

- [10. Spectral Methods, FEM, and MUSCL Schemes](numerics_and_solvers.md)
  - [Example: Spectral 3D Poisson](examples/spectral_3d_poisson.md)
  - [Example: FEM 3D Stokes](examples/fem_3d_stokes.md)
  - [Example: MUSCL Schemes Demo](examples/muscl_schemes_demo.md)
  - [Example: Spectral Performance](examples/spectral_performance.md)
- [11. Matrix-Free Operators and Krylov Solvers](matrix_free_operators.md)
  - [Example: Matrix Free Demo](examples/matrix_free_demo.md)
  - [Example: 2D Heat Diffusion](examples/2d_heat_diffusion.md)

# Part V — Geometry, Meshing, and Flow Network Design

- [12. Constructive Solid Geometry for CFD](geometry_and_meshing.md)
  - [Example: CSG Primitives Demo](examples/csg_primitives_demo.md)
  - [Example: CSG Operations](examples/csg_operations.md)
  - [Example: CSG CFD Simulation](examples/csg_cfd_simulation.md)
  - [Example: Mesh 3D Integration](examples/mesh_3d_integration.md)
- [13. Microfluidic Schematics and Channel Networks](crate_schematics.md)
  - [Example: Bifurcation Schematic](examples/bifurcation_schematic.md)
  - [Example: Venturi Schematic](examples/venturi_schematic.md)
  - [Example: Serpentine Mixing Schematic](examples/serpentine_mixing_schematic.md)
  - [Example: Dimension Scenarios Plots](examples/dimension_scenarios_plots.md)
  - [Example: Comprehensive Arc Demo](examples/comprehensive_arc_demo.md)
  - [Example: Serpentine Venturi Demo](examples/serpentine_venturi_demo.md)
  - [Example: Frustum Channel Demo](examples/frustum_channel_demo.md)
  - [Example: Shell Cuboid Demo](examples/shell_cuboid_demo.md)

# Part VI — Three-Dimensional Flows

- [14. 3-D Navier-Stokes: Bifurcations, Cavitation, and Dean Vortices](crate_3d_flows.md)
  - [Example: Bifurcation 3D Blood](examples/bifurcation_3d_blood.md)
  - [Example: Venturi 3D Cavitation](examples/venturi_3d_cavitation.md)
  - [Example: Serpentine 3D Dean](examples/serpentine_3d_dean.md)
  - [Example: Spectral Poisson 3D](examples/spectral_poisson_3d_crate.md)
- [15. 2-D Flows and Schematic Integration](schematic_integration_2d.md)
  - [Example: Blood Venturi](examples/blood_venturi.md)
  - [Example: TPMS Blood 2D](examples/tpms_blood_2d.md)
  - [Example: Serpentine Venturi 1D vs 2D](examples/serpentine_venturi_1d_vs_2d.md)
  - [Example: Schematic Demo Integration](examples/schematic_demo_integration.md)
  - [Example: Geometry Integration Demo](examples/geometry_integration_demo.md)

# Part VII — Validation and Benchmarking

- [16. Canonical Flow Benchmarks](validation.md)
  - [Example: Cavity Validation](examples/cavity_validation.md)
  - [Example: Pipe Flow Validation](examples/pipe_flow_validation.md)
  - [Example: Comprehensive Validation Suite](examples/comprehensive_validation_suite.md)
  - [Example: Richardson Convergence](examples/richardson_convergence.md)
  - [Example: Blood Poiseuille 2D](examples/blood_poiseuille_2d.md)
- [17. GPU Detection and Performance Profiling](performance_and_atlas.md)
  - [Example: GPU Detection](examples/gpu_detection.md)
  - [Example: SIMD Performance Benchmark](examples/simd_performance_benchmark.md)
  - [Example: Venturi Validated](examples/venturi_validated.md)
- [18. Validation Suite](crate_validation.md)

# Part VIII — Optimization and Device Design

- [19. Multi-Objective Optimization](crate_optim.md)
  - [Example: Cell Sep Audit](examples/cell_sep_audit.md)
  - [Example: Milestone 12 Validation](examples/milestone12_validation.md)
  - [Example: Milestone 12 Option 1](examples/milestone12_option1.md)
  - [Example: Milestone 12 Option 2](examples/milestone12_option2.md)
  - [Example: Milestone 12 GA](examples/milestone12_ga.md)

# Part IX — Atlas Stack Migration

- [20. Migration Overview](migration_overview.md)
- [21. Eunomia: Numeric Traits](migration_eunomia.md)
- [22. Leto: Arrays and Linalg](migration_arrays.md)
- [23. Leto: Geometry](migration_geometry.md)
- [24. Hermes: SIMD Lanes](migration_simd.md)
- [25. Mnemosyne and Themis: Memory](migration_memory.md)
- [26. Moirai: Concurrency](migration_concurrency.md)
- [27. Apollo: FFT](migration_fft.md)
- [28. Leto: GAT Tiling](migration_gat_tiles.md)
- [29. Migration Validation](migration_validation.md)

# Appendix

- [A. Atlas Crate Dependency Map](appendix_dependencies.md)
- [B. Glossary](appendix_glossary.md)
- [C. Changelog](appendix_changelog.md)
- [D. Atlas Migration Reference](appendix_migration.md)
- [E. Book Organization](BOOK_ORGANIZATION.md)
- [F. Parity Artefacts Archive](parity_archive.md)
