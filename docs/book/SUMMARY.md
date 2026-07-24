# Summary

[Introduction](README.md)

---

# Part I — Governing Equations and Physical Models

- [1. CFDrs Architecture and Problem Setup](foundations.md)
  - [Example: cfd_demo](examples/cfd_demo.md)
  - [Example: enhanced_cfd_demo](examples/enhanced_cfd_demo.md)
- [2. Governing Equations and Discretization](governing_equations.md)
- [3. Pressure-Velocity Coupling and Time Integration](pressure_velocity.md)
  - [Example: simplec_pimple_demo](examples/simplec_pimple_demo.md)
  - [Example: adaptive_time_stepping_demo](examples/adaptive_time_stepping_demo.md)

---

# Part II — Turbulence and Multiphase Physics

- [4. Turbulence Models: k-ε, k-ω SST, Smagorinsky LES](turbulence_multiphase.md)
  - [Example: turbulence_models_demo](examples/turbulence_models_demo.md)
  - [Example: turbulence_validation_demo](examples/turbulence_validation_demo.md)
  - [Example: turbulence_momentum_integration_demo](examples/turbulence_momentum_integration_demo.md)
  - [Example: turbulent_channel_flow](examples/turbulent_channel_flow.md)
  - [Figure: Semi-log residual convergence (SIMPLE vs PISO)](figures/residual_convergence_semilog.svg)
- [5. Cavitation: Liquid-Vapour Phase Transition](cavitation.md)
  - [Example: venturi_cavitation](examples/venturi_cavitation.md)
  - [Example: simple_cavitation](examples/simple_cavitation.md)
  - [Example: cavitation_damage_simulation](examples/cavitation_damage_simulation.md)
  - [Example: 1d_venturi_cavitation](examples/1d_venturi_cavitation.md)

---

# Part III — Biomedical and Rheological Flows

- [6. Non-Newtonian Fluids and Blood Rheology](biomedical_flows.md)
  - [Example: blood_rheology_models](examples/blood_rheology_models.md)
  - [Example: blood_flow_1d_validation](examples/blood_flow_1d_validation.md)
- [7. Vascular Bifurcations and Stenosis](vascular_bifurcations.md)
  - [Example: bifurcation_2d_blood_validation](examples/bifurcation_2d_blood_validation.md)
  - [Example: venturi_blood_flow_validation](examples/venturi_blood_flow_validation.md)
  - [Example: cross_fidelity_branching](examples/cross_fidelity_branching.md)
- [8. Microfluidics and Millifluidic Networks](crate_1d_flows.md)
  - [Example: blood_bifurcation](examples/blood_bifurcation.md)
  - [Example: microfluidic_chip](examples/microfluidic_chip.md)
  - [Example: fda_shear_limit_screening](examples/fda_shear_limit_screening.md)
  - [Example: tpms_blood_1d](examples/tpms_blood_1d.md)
  - [Example: medical_millifluidic_screening](examples/medical_millifluidic_screening.md)
  - [Example: hemolysis_serpentine_analysis](examples/hemolysis_serpentine_analysis.md)
  - [Example: serpentine_mixing_comprehensive](examples/serpentine_mixing_comprehensive.md)

---

# Part IV — Numerical Methods and Linear Solvers

- [9. Spectral Methods, FEM, and MUSCL Schemes](numerics_and_solvers.md)
  - [Example: spectral_3d_poisson](examples/spectral_3d_poisson.md)
  - [Example: fem_3d_stokes](examples/fem_3d_stokes.md)
  - [Example: muscl_schemes_demo](examples/muscl_schemes_demo.md)
  - [Example: spectral_performance](examples/spectral_performance.md)
- [10. Matrix-Free Operators and Krylov Solvers](matrix_free_operators.md)
  - [Example: matrix_free_demo](examples/matrix_free_demo.md)
  - [Example: 2d_heat_diffusion](examples/2d_heat_diffusion.md)

---

# Part V — Geometry, Meshing, and Flow Network Design

- [11. Constructive Solid Geometry for CFD](geometry_and_meshing.md)
  - [Example: csg_primitives_demo](examples/csg_primitives_demo.md)
  - [Example: csg_operations](examples/csg_operations.md)
  - [Example: csg_cfd_simulation](examples/csg_cfd_simulation.md)
  - [Example: mesh_3d_integration](examples/mesh_3d_integration.md)
  - [Figure: Structured 16 x 4 channel mesh layout](figures/channel_mesh_layout.svg)
- [12. Microfluidic Schematics and Channel Networks](crate_schematics.md)
  - [Example: bifurcation_schematic](examples/bifurcation_schematic.md)
  - [Example: venturi_schematic](examples/venturi_schematic.md)
  - [Example: serpentine_mixing_schematic](examples/serpentine_mixing_schematic.md)
  - [Example: dimension_scenarios_plots](examples/dimension_scenarios_plots.md)
  - [Example: comprehensive_arc_demo](examples/comprehensive_arc_demo.md)
  - [Example: serpentine_venturi_demo](examples/serpentine_venturi_demo.md)
  - [Example: frustum_channel_demo](examples/frustum_channel_demo.md)
  - [Example: shell_cuboid_demo](examples/shell_cuboid_demo.md)

---

# Part VI — Three-Dimensional Flows

- [13. 3-D Navier-Stokes: Bifurcations, Cavitation, and Dean Vortices](crate_3d_flows.md)
  - [Example: bifurcation_3d_blood](examples/bifurcation_3d_blood.md)
  - [Example: venturi_3d_cavitation](examples/venturi_3d_cavitation.md)
  - [Example: serpentine_3d_dean](examples/serpentine_3d_dean.md)
  - [Example: spectral_poisson_3d (cfd-3d)](examples/spectral_poisson_3d_crate.md)
- [14. 2-D Flows and Schematic Integration](schematic_integration_2d.md)
  - [Example: blood_venturi](examples/blood_venturi.md)
  - [Example: tpms_blood_2d](examples/tpms_blood_2d.md)
  - [Example: serpentine_venturi_1d_vs_2d](examples/serpentine_venturi_1d_vs_2d.md)
  - [Example: schematic_demo_integration](examples/schematic_demo_integration.md)
  - [Example: geometry_integration_demo](examples/geometry_integration_demo.md)

---

# Part VII — Validation and Benchmarking

- [15. Canonical Flow Benchmarks](validation.md)
  - [Example: cavity_validation](examples/cavity_validation.md)
  - [Figure: Lid-driven cavity streamfunction contours (Re ~ 100)](figures/cavity_streamfunction_contour.svg)
  - [Example: pipe_flow_validation](examples/pipe_flow_validation.md)
  - [Figure: Plane Poiseuille parabolic profile (1-D analytical)](figures/poiseuille_parabolic_profile.svg)
  - [Figure: Pipe-flow Reynolds regime map](figures/reynolds_regime_map.svg)
  - [Example: comprehensive_validation_suite](examples/comprehensive_validation_suite.md)
  - [Example: richardson_convergence](examples/richardson_convergence.md)
  - [Figure: Richardson extrapolation log-log (2nd-order reference)](figures/richardson_loglog.svg)
  - [Example: blood_poiseuille_2d](examples/blood_poiseuille_2d.md)
- [16. GPU Detection and Performance Profiling](performance_and_atlas.md)
  - [Example: gpu_detection](examples/gpu_detection.md)
  - [Example: simd_performance_benchmark](examples/simd_performance_benchmark.md)
  - [Example: venturi_validated](examples/venturi_validated.md)
  - [Figure: CFDrs layered architecture on the Atlas stack](figures/architecture_stack.svg)

---

# Part VIII — Optimization and Device Design

- [17. Multi-Objective Optimization](crate_optim.md)
  - [Example: cell_sep_audit](examples/cell_sep_audit.md)
  - [Example: milestone12_validation](examples/milestone12_validation.md)
  - [Example: milestone12_option1](examples/milestone12_option1.md)
  - [Example: milestone12_option2](examples/milestone12_option2.md)
  - [Example: milestone12_ga](examples/milestone12_ga.md)

---

# Appendix

- [A. Atlas Crate Dependency Map](appendix_dependencies.md)
- [B. Atlas Stack Reference: ndarray/nalgebra → leto/hephaestus/coeus](appendix_migration.md)
- [C. Glossary](appendix_glossary.md)
- [D. Changelog](appendix_changelog.md)
- [E. Book Organization](BOOK_ORGANIZATION.md)
- [F. Parity Artefacts Archive (CI Gate Evidence)](../../../parity_artefacts/INDEX.md)

- [Stray test figure (drift fixture)](figures/atlas_drift_fixture_test_only_no_such_figure.svg)
