//! 1D Blood Flow Validation Tests
//!
//! Validation of 1D network solver for non-Newtonian blood flow in:
//! 1. Stenosis (Venturi-like varying cross-section)
//! 2. Bifurcation (Y-junction)
//!
//! References:
//! - Gijsen, F. J. H., et al. (1999). "A new experimental technique to distinguish between yield stress and shear-thinning effects of blood."
//!   Journal of Biomechanics, 32(6), 617-623.
//! - Ku, D. N. (1997). "Blood flow in arteries." Annual review of fluid mechanics, 29(1), 399-434.

use super::{LiteratureValidation, ValidationReport};
use crate::scalar::{self, ValidationScalar};
use cfd_1d::{
    // Channel geometry types — re-exported at crate root from domain::channel
    ChannelGeometry,
    ChannelType,
    // Network types — re-exported at crate root from domain::network
    ComponentType,
    CrossSection,
    EdgeProperties,
    Network,
    NetworkBuilder,
    // Solver types — re-exported at crate root from solver::core
    NetworkProblem,
    NetworkSolver,
    SolverConfig,
    SurfaceProperties,
    Wettability,
};
use cfd_core::{error::Result, physics::fluid::non_newtonian::CarreauYasuda};
use std::collections::HashMap;

/// Stenosis validation case (Venturi-like geometry)
///
/// Simulates flow through a channel with a constriction:
/// Inlet -> Wide (D=10mm) -> Narrow (D=5mm) -> Wide (D=10mm) -> Outlet
///
/// Validates:
/// - Mass conservation
/// - Pressure drop increase in constriction
/// - Shear-thinning behavior (lower viscosity in constriction)
pub struct StenosisValidation<T: ValidationScalar> {
    accuracy: T,
}

impl<T: ValidationScalar> StenosisValidation<T> {
    /// Create a new stenosis validation test
    pub fn new() -> Self {
        Self {
            accuracy: scalar::from_f64::<T>(1e-2), // 1% accuracy target
        }
    }
}

impl<T: ValidationScalar> LiteratureValidation<T> for StenosisValidation<T> {
    fn validate(&self) -> Result<ValidationReport<T>> {
        let mut builder = NetworkBuilder::new();

        // 1. Define nodes
        let n_inlet = builder.add_inlet("Inlet".to_string());
        let n1 = builder.add_junction("J1".to_string());
        let n2 = builder.add_junction("J2".to_string());
        let n3 = builder.add_junction("J3".to_string());
        let n_outlet = builder.add_outlet("Outlet".to_string());

        // 2. Define edges (channels)
        // Segment 1: Wide (D=10mm, L=50mm)
        let e1 = builder.connect_with_pipe(n_inlet, n1, "Seg1_Wide".to_string());
        // Segment 2: Converging (modeled as intermediate pipe, D=7.5mm, L=10mm) - Simplified 1D
        let e2 = builder.connect_with_pipe(n1, n2, "Seg2_Converge".to_string());
        // Segment 3: Narrow (Stenosis, D=5mm, L=20mm)
        let e3 = builder.connect_with_pipe(n2, n3, "Seg3_Stenosis".to_string());
        // Segment 4: Diverging (modeled as intermediate pipe, D=7.5mm, L=10mm)
        let e4 = builder.connect_with_pipe(n3, n_outlet, "Seg4_Diverge".to_string());

        // 3. Build Graph
        let graph = builder.build()?;

        // 4. Create Network with Blood
        // Carreau-Yasuda model for blood
        let blood = CarreauYasuda::<T>::blood();
        let mut network = Network::new(graph, blood);

        // 5. Assign Properties (Geometry)
        let mm = scalar::from_f64::<T>(1e-3);
        let d_wide = scalar::from_f64::<T>(10.0) * mm;
        let d_mid = scalar::from_f64::<T>(7.5) * mm;
        let d_narrow = scalar::from_f64::<T>(5.0) * mm;

        // Ensure L/D > 10 for Hagen-Poiseuille validity
        let l_wide = scalar::from_f64::<T>(150.0) * mm;
        let l_trans = scalar::from_f64::<T>(100.0) * mm;
        let l_narrow = scalar::from_f64::<T>(60.0) * mm;

        // Helper to create properties
        let create_props = |d: T, l: T| -> EdgeProperties<T> {
            let pi = scalar::from_f64::<T>(std::f64::consts::PI);
            let area = pi * d * d / scalar::from_f64::<T>(4.0);
            // Estimate initial resistance using Poiseuille for water (approx) to avoid huge initial flows
            // R = 128 * mu * L / (pi * D^4). mu_water ~ 1e-3.
            let mu = scalar::from_f64::<T>(1e-3);
            let r_init = scalar::from_f64::<T>(128.0) * mu * l / (pi * scalar::powi(d, 4));

            EdgeProperties {
                id: String::new(), // Overwritten by network map
                resistance_update_policy: cfd_1d::ResistanceUpdatePolicy::FlowDependent,
                component_type: ComponentType::Pipe,
                length: l,
                area,
                hydraulic_diameter: Some(d),
                resistance: r_init,
                geometry: Some(ChannelGeometry {
                    channel_type: ChannelType::Straight,
                    length: l,
                    cross_section: CrossSection::Circular { diameter: d },
                    surface: SurfaceProperties {
                        roughness: T::zero(),
                        contact_angle: None,
                        surface_energy: None,
                        wettability: Wettability::Hydrophilic,
                    },
                    variations: Vec::new(),
                }),
                properties: HashMap::new(),
            }
        };

        network.add_edge_properties(e1, create_props(d_wide, l_wide));
        network.add_edge_properties(e2, create_props(d_mid, l_trans));
        network.add_edge_properties(e3, create_props(d_narrow, l_narrow));
        network.add_edge_properties(e4, create_props(d_mid, l_trans));

        // Ensure resistances are initialized from properties before solving
        network.update_resistances()?;

        // 6. Set Boundary Conditions
        // Reduce pressure to ensure Laminar flow for shear thinning validation
        let p_in = scalar::from_f64::<T>(100.0);
        let p_out = T::zero();
        network.set_pressure(n_inlet, p_in);
        network.set_pressure(n_outlet, p_out);

        // 7. Solve
        let problem = NetworkProblem::new(network);
        let solver = NetworkSolver::with_config(SolverConfig {
            tolerance: scalar::from_f64::<T>(1e-6),
            max_iterations: 100,
        });

        let solution = solver.solve_network(&problem)?;

        // 8. Analysis and Validation
        let q1 = scalar::abs(
            solution
                .flow_rates
                .get(e1.index())
                .copied()
                .unwrap_or(T::zero()),
        );
        let q3 = scalar::abs(
            solution
                .flow_rates
                .get(e3.index())
                .copied()
                .unwrap_or(T::zero()),
        );
        let q_err = scalar::abs(q1 - q3);

        let dp1 = scalar::abs(solution.pressures[n_inlet.index()] - solution.pressures[n1.index()]);
        let dp3 = scalar::abs(solution.pressures[n2.index()] - solution.pressures[n3.index()]);

        let grad1 = dp1 / l_wide;
        let grad3 = dp3 / l_narrow;

        let grad_ratio = grad3 / grad1;

        let mu_ratio = grad_ratio / scalar::from_f64::<T>(16.0);
        // Check if apparent viscosity is significantly reduced compared to Newtonian baseline
        // For Carreau-Yasuda blood model, high shear in stenosis drastically lowers viscosity.
        let is_viscosity_reduced = mu_ratio < scalar::from_f64::<T>(0.95);

        // Debug info from edges
        let mut details = String::new();
        use petgraph::visit::EdgeRef;
        for edge_ref in solution.graph.edge_references() {
            let w = edge_ref.weight();
            details.push_str(&format!(
                "Edge {} R: {:e}, k: {:e}\n",
                edge_ref.id().index(),
                scalar::to_f64(w.resistance),
                scalar::to_f64(w.quad_coeff)
            ));
        }

        details.push_str(&format!("Flow Rate Q: {:e} m3/s\n", scalar::to_f64(q1)));
        details.push_str(&format!(
            "Pressure Gradient Wide: {:e} Pa/m\n",
            scalar::to_f64(grad1)
        ));
        details.push_str(&format!(
            "Pressure Gradient Stenosis: {:e} Pa/m\n",
            scalar::to_f64(grad3)
        ));
        details.push_str(&format!(
            "Gradient Ratio: {:.2} (Newtonian exp: 16.0)\n",
            scalar::to_f64(grad_ratio)
        ));
        details.push_str(&format!(
            "Apparent Viscosity Ratio (Stenosis/Wide): {:.2}\n",
            scalar::to_f64(mu_ratio)
        ));

        // In turbulent regime (if k > 0), simple ratio checks fail. We accept if Q is reasonable and mass conserved.
        // Or if we specifically targeted laminar regime (which we did by lowering pressure), we expect shear thinning.
        // However, if transition still happens or if "shear thinning" effect is small, we should be lenient.
        // We mainly validate that the simulation runs physically correctly AND shows shear thinning.
        let passed = q_err < scalar::from_f64::<T>(1e-10) && is_viscosity_reduced;

        Ok(ValidationReport {
            test_name: "1D Stenosis Blood Flow".to_string(),
            citation: self.citation().to_string(),
            max_error: q_err,
            avg_error: q_err,
            passed,
            details,
        })
    }

    fn citation(&self) -> &'static str {
        "Ku, D. N. (1997). Blood flow in arteries."
    }

    fn expected_accuracy(&self) -> T {
        self.accuracy
    }
}

/// Bifurcation validation case (Y-Junction)
///
/// Simulates flow splitting into two branches.
/// Inlet -> Parent -> Junction -> (Branch 1, Branch 2)
///
/// Validates:
/// - Flow splitting based on resistance
/// - Mass conservation at junction
pub struct BifurcationValidation<T: ValidationScalar> {
    accuracy: T,
}

impl<T: ValidationScalar> BifurcationValidation<T> {
    /// Create a new bifurcation validation test
    pub fn new() -> Self {
        Self {
            accuracy: scalar::from_f64::<T>(1e-2),
        }
    }
}

impl<T: ValidationScalar> LiteratureValidation<T> for BifurcationValidation<T> {
    fn validate(&self) -> Result<ValidationReport<T>> {
        let mut builder = NetworkBuilder::new();

        let n_inlet = builder.add_inlet("Inlet".to_string());
        let n_junc = builder.add_junction("Bifurcation".to_string());
        let n_out1 = builder.add_outlet("Out1".to_string());
        let n_out2 = builder.add_outlet("Out2".to_string());

        let e_parent = builder.connect_with_pipe(n_inlet, n_junc, "Parent".to_string());
        let e_branch1 = builder.connect_with_pipe(n_junc, n_out1, "Branch1".to_string());
        let e_branch2 = builder.connect_with_pipe(n_junc, n_out2, "Branch2".to_string());

        let graph = builder.build()?;
        let blood = CarreauYasuda::<T>::blood();
        let mut network = Network::new(graph, blood);

        let mm = scalar::from_f64::<T>(1e-3);
        // Parent: 10mm dia, 150mm len (L/D > 10)
        // Branch 1: 8mm dia, 100mm len
        // Branch 2: 6mm dia, 100mm len (Higher resistance)

        let create_props = |d: T, l: T| -> EdgeProperties<T> {
            let pi = scalar::from_f64::<T>(std::f64::consts::PI);
            let area = pi * d * d / scalar::from_f64::<T>(4.0);
            // Initial resistance guess
            let mu = scalar::from_f64::<T>(1e-3);
            let r_init = scalar::from_f64::<T>(128.0) * mu * l / (pi * scalar::powi(d, 4));

            EdgeProperties {
                id: String::new(),
                resistance_update_policy: cfd_1d::ResistanceUpdatePolicy::FlowDependent,
                component_type: ComponentType::Pipe,
                length: l,
                area,
                hydraulic_diameter: Some(d),
                resistance: r_init,
                geometry: Some(ChannelGeometry {
                    channel_type: ChannelType::Straight,
                    length: l,
                    cross_section: CrossSection::Circular { diameter: d },
                    surface: SurfaceProperties {
                        roughness: T::zero(),
                        contact_angle: None,
                        surface_energy: None,
                        wettability: Wettability::Hydrophilic,
                    },
                    variations: Vec::new(),
                }),
                properties: HashMap::new(),
            }
        };

        network.add_edge_properties(
            e_parent,
            create_props(
                scalar::from_f64::<T>(10.0) * mm,
                scalar::from_f64::<T>(150.0) * mm,
            ),
        );
        network.add_edge_properties(
            e_branch1,
            create_props(
                scalar::from_f64::<T>(8.0) * mm,
                scalar::from_f64::<T>(100.0) * mm,
            ),
        );
        network.add_edge_properties(
            e_branch2,
            create_props(
                scalar::from_f64::<T>(6.0) * mm,
                scalar::from_f64::<T>(100.0) * mm,
            ),
        );

        // Ensure resistances are initialized
        network.update_resistances()?;

        // Low pressure for laminar flow
        network.set_pressure(n_inlet, scalar::from_f64::<T>(10.0));
        network.set_pressure(n_out1, T::zero());
        network.set_pressure(n_out2, T::zero());

        let problem = NetworkProblem::new(network);
        let solver = NetworkSolver::with_config(SolverConfig {
            tolerance: scalar::from_f64::<T>(1e-6),
            max_iterations: 100,
        });

        let solution = solver.solve_network(&problem)?;

        let q_parent = scalar::abs(
            solution
                .flow_rates
                .get(e_parent.index())
                .copied()
                .unwrap_or(T::zero()),
        );
        let q_b1 = scalar::abs(
            solution
                .flow_rates
                .get(e_branch1.index())
                .copied()
                .unwrap_or(T::zero()),
        );
        let q_b2 = scalar::abs(
            solution
                .flow_rates
                .get(e_branch2.index())
                .copied()
                .unwrap_or(T::zero()),
        );

        // 1. Mass Conservation
        let sum_out = q_b1 + q_b2;
        let mass_err = scalar::abs(q_parent - sum_out);

        // 2. Flow Splitting
        let flow_ratio = q_b1 / q_b2;
        let expected_ratio_newtonian = scalar::from_f64::<T>(3.16);

        // Shear rate analysis:
        // R ~ 1/D^4. Q ~ D^4. V ~ Q/D^2 ~ D^2. Shear ~ V/D ~ D.
        // So in parallel branches, the WIDER branch has HIGHER shear rate!
        // Higher shear -> Lower viscosity -> Lower resistance -> Higher flow.
        // So Q1 (wide) increases relative to Q2 (narrow).
        // Expect flow_ratio > 3.16.
        let is_shear_thinning_split = flow_ratio > expected_ratio_newtonian;

        let mut details = String::new();
        details.push_str(&format!("Q Parent: {:e}\n", scalar::to_f64(q_parent)));
        details.push_str(&format!("Q Branch 1 (8mm): {:e}\n", scalar::to_f64(q_b1)));
        details.push_str(&format!("Q Branch 2 (6mm): {:e}\n", scalar::to_f64(q_b2)));
        details.push_str(&format!(
            "Flow Ratio Q1/Q2: {:.2} (Newtonian exp: 3.16)\n",
            scalar::to_f64(flow_ratio)
        ));
        details.push_str(&format!("Mass Error: {:e}\n", scalar::to_f64(mass_err)));

        let passed = mass_err < scalar::from_f64::<T>(1e-10) && is_shear_thinning_split;

        Ok(ValidationReport {
            test_name: "1D Bifurcation Blood Flow".to_string(),
            citation: "Standard Fluid Dynamics".to_string(),
            max_error: mass_err,
            avg_error: mass_err,
            passed,
            details,
        })
    }

    fn citation(&self) -> &'static str {
        "Standard Fluid Dynamics - Bifurcation"
    }

    fn expected_accuracy(&self) -> T {
        self.accuracy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stenosis_blood_flow() {
        let validation = StenosisValidation::<f64>::new();
        let report = validation.validate().expect("Validation failed to run");
        println!("{}", report.details);
        assert!(
            report.passed,
            "Stenosis validation failed: {}",
            report.details
        );
    }

    #[test]
    fn test_bifurcation_blood_flow() {
        let validation = BifurcationValidation::<f64>::new();
        let report = validation.validate().expect("Validation failed to run");
        println!("{}", report.details);
        assert!(
            report.passed,
            "Bifurcation validation failed: {}",
            report.details
        );
    }
}
