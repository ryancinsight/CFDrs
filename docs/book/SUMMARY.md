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
  - [Example: turbulence_momentum_integration_demo](examples/turbulence_momentum_integration_demo.md)
  - [Example: venturi_cavitation](examples/venturi_cavitation.md)
  - [Example: simple_cavitation](examples/simple_cavitation.md)
  - [Example: cavitation_damage_simulation](examples/cavitation_damage_simulation.md)
  - [Example: 1d_venturi_cavitation](examples/1d_venturi_cavitation.md)

# Part IV — Biomedical and Specialized Flows

- [Blood Flow and Rheology Workflows](biomedical_flows.md)
  - [Example: blood_flow_1d_validation](examples/blood_flow_1d_validation.md)
  - [Example: bifurcation_2d_blood_validation](examples/bifurcation_2d_blood_validation.md)
  - [Example: blood_rheology_models](examples/blood_rheology_models.md)
  - [Example: venturi_blood_flow_validation](examples/venturi_blood_flow_validation.md)
  - [Example: serpentine_mixing_comprehensive](examples/serpentine_mixing_comprehensive.md)
  - [Example: cross_fidelity_branching](examples/cross_fidelity_branching.md)

# Part V — Discretization and Solvers

- [Spectral, FEM, MUSCL, and Matrix-Free Methods](numerics_and_solvers.md)
  - [Example: spectral_3d_poisson](examples/spectral_3d_poisson.md)
  - [Example: fem_3d_stokes](examples/fem_3d_stokes.md)
  - [Example: matrix_free_demo](examples/matrix_free_demo.md)
  - [Example: muscl_schemes_demo](examples/muscl_schemes_demo.md)
  - [Example: simplec_pimple_demo](examples/simplec_pimple_demo.md)
  - [Example: 2d_heat_diffusion](examples/2d_heat_diffusion.md)

# Part VI — Geometry, Meshing, and CSG

- [Geometry Construction and CFD Coupling](geometry_and_meshing.md)
  - [Example: csg_primitives_demo](examples/csg_primitives_demo.md)
  - [Example: csg_operations](examples/csg_operations.md)
  - [Example: csg_cfd_simulation](examples/csg_cfd_simulation.md)
  - [Example: mesh_3d_integration](examples/mesh_3d_integration.md)
  - [Example: dimension_scenarios_plots](examples/dimension_scenarios_plots.md)

# Part VII — Performance and Atlas Integration

- [SIMD, GPU, and Backend Migration](performance_and_atlas.md)
  - [Example: simd_performance_benchmark](examples/simd_performance_benchmark.md)
  - [Example: gpu_detection](examples/gpu_detection.md)
  - [Example: spectral_performance](examples/spectral_performance.md)
  - [Example: venturi_validated](examples/venturi_validated.md)

---

# Part VIII — Crate-Level Examples

- [Validation Suite](crate_validation.md)
  - [Example: comprehensive_validation_suite](examples/comprehensive_validation_suite.md)
  - [Example: richardson_convergence](examples/richardson_convergence.md)
  - [Example: blood_poiseuille_2d](examples/blood_poiseuille_2d.md)

- [3-D Flows](crate_3d_flows.md)
  - [Example: spectral_poisson_3d (cfd-3d)](examples/spectral_poisson_3d_crate.md)
  - [Example: bifurcation_3d_blood](examples/bifurcation_3d_blood.md)
  - [Example: venturi_3d_cavitation](examples/venturi_3d_cavitation.md)
  - [Example: serpentine_3d_dean](examples/serpentine_3d_dean.md)

- [1-D Biomedical Flows](crate_1d_flows.md)
  - [Example: blood_bifurcation (cfd-1d)](examples/blood_bifurcation.md)
  - [Example: microfluidic_chip](examples/microfluidic_chip.md)
  - [Example: fda_shear_limit_screening](examples/fda_shear_limit_screening.md)
  - [Example: venturi_parallel_analysis](examples/venturi_parallel_analysis.md)
  - [Example: tpms_blood_1d](examples/tpms_blood_1d.md)
  - [Example: cavitation_venturi_analysis](examples/cavitation_venturi_analysis.md)
  - [Example: medical_millifluidic_screening](examples/medical_millifluidic_screening.md)
  - [Example: hemolysis_serpentine_analysis](examples/hemolysis_serpentine_analysis.md)

- [2-D and Schematic Examples](crate_schematics.md)
  - [Example: bifurcation_schematic](examples/bifurcation_schematic.md)
  - [Example: venturi_schematic](examples/venturi_schematic.md)
  - [Example: serpentine_mixing_schematic](examples/serpentine_mixing_schematic.md)
  - [Example: blood_venturi](examples/blood_venturi.md)
  - [Example: tpms_blood_2d](examples/tpms_blood_2d.md)
  - [Example: serpentine_venturi_1d_vs_2d](examples/serpentine_venturi_1d_vs_2d.md)
  - [Example: schematic_demo_integration](examples/schematic_demo_integration.md)
  - [Example: geometry_integration_demo](examples/geometry_integration_demo.md)

- [Optimization](crate_optim.md)
  - [Example: cell_sep_audit](examples/cell_sep_audit.md)
  - [Example: milestone12_validation](examples/milestone12_validation.md)
  - [Example: milestone12_report](examples/milestone12_report.md)
  - [Example: milestone12_option1](examples/milestone12_option1.md)
  - [Example: milestone12_option2](examples/milestone12_option2.md)
  - [Example: milestone12_ga](examples/milestone12_ga.md)

---

# Appendix

- [Atlas Crate Dependency Map](appendix_dependencies.md)
- [Migration Notes: ndarray/nalgebra/burn → leto/hephaestus/coeus](appendix_migration.md)
